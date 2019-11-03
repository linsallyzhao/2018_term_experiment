#include <chrono>
#include <math.h>
#include <iostream>
#include <random>
#include <fstream>
struct Option {
    double inRate, pricingVol, realVol, expire, stock, strike;
    int steps;
    double dt, discount_factor;
};

double get_dOne(const Option &option, const bool &real, const int &stepNow);
double one_step(const double &old, const Option &option);


int main(void) {
    Option option;
    option.inRate = 0.02;
    option.pricingVol = 0.1;
    option.realVol = 0.3;
    option.stock = 100; //beginning price of stock
    option.strike = 100;
    option.steps = 100;//how many hedging steps we will record
    option.expire = 1.0;
    option.dt = option.expire / option.steps;
    //This is ONLY suitable for the expiry discount
    option.discount_factor = exp(-option.inRate * option.expire);
    int npaths = 1000000;

    double dOne_1 = get_dOne(option, true, 0); //d_1 calculated with real vol
    double dOne_2 = get_dOne(option, false, 0); //d_1 calculate with pricing vol
    double NdOne_1 = std::erfc(-dOne_1/std::sqrt(2))/2;//normal CDF dOne_1
    double NdOne_2 = std::erfc(-dOne_2/std::sqrt(2))/2;
    double initialDelta_1 = option.discount_factor*NdOne_1;
    //Normal CDF dOne_1 times discount factor
    double initialDelta_2 = option.discount_factor*NdOne_2;
    //Normal CDF dOne_2 times discount factor
    double initialOptionPrice1 = option.stock*NdOne_1 - option.strike*option.discount_factor
        *std::erfc(-(dOne_1-(option.realVol*pow(option.expire, 0.5)))/std::sqrt(2))/2;
    double initialOptionPrice2 = option.stock*NdOne_2 - option.strike*option.discount_factor
        *std::erfc(-(dOne_2-(option.pricingVol*pow(option.expire, 0.5)))/std::sqrt(2))/2;

    std::cout << "price diff of two vols: " << initialOptionPrice1 - initialOptionPrice2 << std::endl;

    double currentHedge_1, currentHedge_2, cash_1, cash_2, oldPrice_1, delta_1, delta_2, payoff,
           oldPrice_2, optionPrice1, optionPrice2, hedgePosition_2, hedgePosition_1, oldStock,
           changeStock, localDF, dailyClose_1, dailyClose_2;
    std::ofstream finalCashEachDay;
    finalCashEachDay.open("FC_100_1MM");
    finalCashEachDay << "finalCash_realVol,finalCash_pricingVol" << std::endl;
    std::ofstream dailyHedgePosition;
    dailyHedgePosition.open("DH_100_1MM");
    std::ofstream dailyCloseCash;
    dailyCloseCash.open("DC_100_1MM");
    for (int n = 1; n <= npaths; n++) {
        option.stock = 100.00;
        delta_1 = initialDelta_1;
        delta_2 = initialDelta_2;
        cash_1 = -initialOptionPrice2+delta_1*option.stock;
        cash_2 = -initialOptionPrice2+delta_2*option.stock;

        if (n == 1)
            std::cout << "initial_cash_1 " << cash_1 << " initial_cash_2 " << cash_2 << std::endl;

        optionPrice1 = initialOptionPrice1;
        optionPrice2 = initialOptionPrice2;
        //std::cout << "hedgePosition: ";
        for (int stepNow = 1; stepNow <= option.steps; stepNow++){
            localDF = exp(-option.inRate * (option.expire-stepNow*option.dt));
            oldStock = option.stock;
            option.stock = one_step(oldStock, option);
            changeStock = option.stock-oldStock;
            currentHedge_1 = delta_1;
            currentHedge_2 = delta_2;
            dOne_1 = get_dOne(option, true, stepNow);
            dOne_2 = get_dOne(option, false, stepNow);
            NdOne_1 = std::erfc(-dOne_1/std::sqrt(2))/2;
            NdOne_2 = std::erfc(-dOne_2/std::sqrt(2))/2;
            delta_1 = localDF*NdOne_1;
            delta_2 = localDF*NdOne_2;
            oldPrice_1 = optionPrice1;
            oldPrice_2 = optionPrice2;
            optionPrice1 = option.stock*NdOne_1 - option.strike*localDF
                *std::erfc(-(dOne_1-(option.realVol*pow(option.expire-stepNow*option.dt, 0.5)))/std::sqrt(2))/2;
            optionPrice2 = option.stock*NdOne_2 - option.strike*localDF
                *std::erfc(-(dOne_2-(option.pricingVol*pow(option.expire-stepNow*option.dt, 0.5)))/std::sqrt(2))/2;
            hedgePosition_1 = (optionPrice1-oldPrice_1)-currentHedge_1*changeStock; //ideally, this should be close to zero
            hedgePosition_2 = (optionPrice2-oldPrice_2)-currentHedge_2*changeStock;
            cash_1 += (delta_1-currentHedge_1)*option.stock;
            cash_2 += (delta_2-currentHedge_2)*option.stock;

            dailyClose_1 = cash_1+optionPrice2-delta_1*option.stock;
            dailyClose_2 = cash_2+optionPrice2-delta_2*option.stock;

            //output to a file
            dailyHedgePosition << hedgePosition_1 << "," << hedgePosition_2 << ",";
            dailyCloseCash << dailyClose_1 << "," << dailyClose_2 << ",";
            //////////////////////////////////////

        }
        //output to a file
        dailyHedgePosition << std::endl;
        dailyCloseCash << std::endl;
        //////////////////////////////////
        payoff = std::max(option.stock-option.strike, 0.0);
        cash_1 += payoff-delta_1*option.stock;
        cash_2 += payoff-delta_2*option.stock;
        //output to a file
        finalCashEachDay << cash_1 << "," << cash_2 << std::endl;
        ////////////////////////////////////
    }

    return 0;
}

//Return d_1 in BSE details 2017 page 25
double get_dOne(const Option &option, const bool &real, const int &stepNow){
    double vol;
    if (real){
        vol = option.realVol;
    }
    else{
        vol = option.pricingVol;
    }
    return (log(option.stock/option.strike)+(option.inRate+0.5*pow(vol, 2))
        *(option.expire-stepNow*option.dt))/(vol*pow(option.expire-stepNow*option.dt, 0.5));
}

double one_step(const double &old, const Option &option) {
    std::normal_distribution<double> norDist(0.0, 1.0);
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    double one_draw = norDist(generator);
    double new_price =
        old * (1 + option.inRate * option.dt + option.realVol * sqrt(option.dt) * one_draw);
    return new_price;
}

