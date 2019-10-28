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
    option.stock = 100;
    option.strike = 100;
    option.steps = 252;
    option.expire = 1.0;
    option.dt = option.expire / option.steps;
    option.discount_factor = exp(-option.inRate * option.expire);
    int npaths = 100000;

    double dOne_1 = get_dOne(option, true, 0);
    double dOne_2 = get_dOne(option, false, 0);
    double initialDelta_1 = std::erfc(-dOne_1/std::sqrt(2))/2;
    double initialDelta_2 = std::erfc(-dOne_2/std::sqrt(2))/2;
    double initialOptionPrice1 = option.stock*initialDelta_1 - option.strike*option.discount_factor
        *std::erfc(-(dOne_1-(option.realVol*pow(option.expire, 0.5)))/std::sqrt(2))/2;
    double initialOptionPrice2 = option.stock*initialDelta_2 - option.strike*option.discount_factor
        *std::erfc(-(dOne_2-(option.pricingVol*pow(option.expire, 0.5)))/std::sqrt(2))/2;

    std::cout << "price diff of two vols: " << initialOptionPrice1 - initialOptionPrice2 << std::endl;

    double currentHedge_1, currentHedge_2, cash_1, cash_2, oldPrice_1, delta_1, delta_2, payoff,
           oldPrice_2, optionPrice1, optionPrice2, hedgePosition_2, hedgePosition_1, oldStock, changeStock,
           closeToday_1, closeToday_2;
    std::ofstream output;
    output.open("finalCash_hedge3TimesPerDay_1mm");
    output << "finalCash_realVol,finalCash_pricingVol" << std::endl;
    //std::ofstream anotherOutPut;
    //anotherOutPut.open("daylyPositionHedged_hedgeEveryTwoDay_1mm");
    std::ofstream closeTest;
    closeTest.open("daylyClose_252Days");
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
            oldStock = option.stock;
            option.stock = one_step(option.stock, option);
            changeStock = option.stock-oldStock;
            currentHedge_1 = delta_1;
            currentHedge_2 = delta_2;
            dOne_1 = get_dOne(option, true, stepNow);
            dOne_2 = get_dOne(option, false, stepNow);
            delta_1 = std::erfc(-dOne_1/std::sqrt(2))/2;
            delta_2 = std::erfc(-dOne_2/std::sqrt(2))/2;
            //oldPrice_1 = optionPrice1;
            //oldPrice_2 = optionPrice2;
            //optionPrice1 = option.stock*delta_1 - option.strike*option.discount_factor
            //    *std::erfc(-(dOne_1-(option.realVol*pow(option.expire-stepNow*option.dt, 0.5)))/std::sqrt(2))/2;
            optionPrice2 = option.stock*delta_1 - option.strike*option.discount_factor
                *std::erfc(-(dOne_1-(option.pricingVol*pow(option.expire-stepNow*option.dt, 0.5)))/std::sqrt(2))/2;
            //hedgePosition_1 = (optionPrice1-oldPrice_1)-currentHedge_1*changeStock; //ideally, this should be close to zero
            //hedgePosition_2 = (optionPrice2-oldPrice_2)-currentHedge_2*changeStock;

            //std::cout << hedgePosition_1 << " " << hedgePosition_2 << "    ";
            //anotherOutPut << hedgePosition_1 << "," << hedgePosition_2 << ",";

            cash_1 += (delta_1-currentHedge_1)*option.stock;
            cash_2 += (delta_2-currentHedge_2)*option.stock;
            closeToday_1 = cash_1+optionPrice2-delta_1*option.stock;
            closeToday_2 = cash_2+optionPrice2-delta_2*option.stock;
            closeTest << closeToday_1 << "," << closeToday_2 << ",";
        }

        //std::cout << std::endl;
        //anotherOutPut << std::endl;
        closeTest << std::endl;
        payoff = std::max(option.stock-option.strike, 0.0);
        cash_1 += payoff-delta_1*option.stock;
        cash_2 += payoff-delta_2*option.stock;
        //std::cout << "final_cash_1 " << cash_1 << " final_cash_2 " << cash_2 << std::endl;
        output << cash_1 << "," << cash_2 << std::endl;
    }

    return 0;
}

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

