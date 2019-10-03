#include <chrono>
#include <math.h>
#include <iostream>
#include <random>
#include <vector>
#include <unordered_map>
#include <string>
struct Option {
    double inRate, pricingVol, realVol, expire, stock, strike;
    int steps;
    double dt, discount_factor;
};

double get_dOne(const Option &option);
double get_price(const Option &option);
double one_step(const double &old, const Option &option);


int main(void) {
    Option option;
    option.inRate = 0.02;
    option.pricingVol = 0.1;
    option.realVol = 0.3;
    option.stock = 100;
    option.strike = 100;
    int hedgingTimes = 252;
    int changesEachDay = 2;
    option.steps = hedgingTimes * changesEachDay;
    option.expire = hedgingTimes / 252.0;
    option.dt = option.expire / option.steps;
    option.discount_factor = exp(-option.inRate * option.expire);
    int npaths = 10000;
    
    double dOne_1 = get_dOne(option, 1, 0);
    double dOne_2 = get_dOne(option, 0, 0);
    double delta_1 = std::erfc(-dOne_1/std::sqrt(2))/2;
    double delta_2 = std::erfc(-dOne_2/std::sqrt(2))/2;
    double optionPrice1 = option.stock*delta_1 - option.strike*option.discount_factor
        *std::erfc(-(dOne_1-(vol*pow(option.expire-stepNow*option.dt, 0.5)))/std::sqrt(2))/2;
    double optionPrice2 = option.stock*delta_2 - option.strike*option.discount_factor
        *std::erfc(-(dOne_2-(vol*pow(option.expire-stepNow*option.dt, 0.5)))/std::sqrt(2))/2;
    double portfolio_1 = optionPrice1 + delta_1 * option.stock;
    double portfolio_2 = optionPrice2 + delta_2 * option.stock;
    double currentHedge_1, currentHedge_2;
    double oldPrice_1, oldPrice_2;
    for (int n = 1; n <= npaths; n++) {
        option.stock = 100.00;
        for (stepNow = 1; stepNow < option.steps; stepNow++){
            if (stepNow % changesEachDay != 0){
                option.stock = one_step(option.stock, option);
                continue;
            }
            option.stock = one_step(option.stock, option);
            currentHedge_1 = delta_1;
            currentHedge_2 = delta_2;
            dOne_1 = get_dOne(option, 1, stepNow);
            dOne_2 = get_dOne(option, 0, stepNow);
            delta_1 = std::erfc(-dOne_1/std::sqrt(2))/2;
            delta_2 = std::erfc(-dOne_2/std::sqrt(2))/2;
            oldPrice_1 = optionPrice1;
            oldPrice_2 = optionPrice2;
            optionPrice1 = option.stock*delta_1 - option.strike*option.discount_factor
                *std::erfc(-(dOne_1-(vol*pow(option.expire-stepNow*option.dt, 0.5)))/std::sqrt(2))/2;
            optionPrice2 = option.stock*delta_1 - option.strike*option.discount_factor
                *std::erfc(-(dOne_1-(vol*pow(option.expire-stepNow*option.dt, 0.5)))/std::sqrt(2))/2;
            portfolio_1 = (optionPrice1-oldPrice1)+(delta_1-currentHedge_1)*option.stock;
            portfolio_2 = (optionPrice2-oldPrice2)+(delta_2-currentHedge_2)*option.stock;
        }
        //record the match (time decay and delta hedging) and separately record
        //the cost of portfolio on the whole way that means record cash
        //position on the way

    }




    return 0;
}

double get_dOne(const Option &option, const bool &real, const int &stepNow){
    if (real){
        double vol = option.realVol;
    }
    else{
        double vol = option.pricingVol;
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

