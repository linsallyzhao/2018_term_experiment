#include <iostream>
#include <random>
#include <cmath>
#include <chrono>

double one_step (const double &, double &, double &, double &);
double end_price (const int &, double &, double &, double &, double &);
double max_price (const int &, double &, double &, double &, double &);
double min_price (const int &, double &, double &, double &, double &);
double get_payoff(int &, int &, double &, double &, double &, double &, double &);
double diff_simu_min_price(const int &steps, double &start, double &inRate, double &vol, double &dt); 


int main() {
//swear I am going to put these into a strucure
    double inRate = 0.1, sigma = 0.3, expire = 0.5;
    double stock = 1, strike = 1;
    int steps = 1;
    double dt = expire / steps;
    double discount_factor = exp(-inRate * expire);
    int npaths = 1000000;
    double running_avg = 0;
    double payoff;
    for (int type = 0; type <= 1; type++) {
        for (int n = 1; n <= npaths; n++) {
            payoff = get_payoff(type, steps, stock, inRate, sigma, dt, strike);
            running_avg = (running_avg * (n - 1) + payoff) / n;
        }

        double price = discount_factor * running_avg;
        if (type == 0 ) {
            std::cout << "Monte carlo price for lookback min_price put is " << price << std::endl;
        }
        else if (type == 1) {
            std::cout << "Monte carlo price for lookback diff_simu put is " << price << std:: endl;
        }
    }

    return 0;
}


double get_payoff(int &type, int &steps, double &stock, double &inRate,
            double &sigma, double &dt, double &strike) {
    double payoff;

    switch(type) {
        case 0 : payoff = std::max(strike - min_price(steps, stock, inRate, sigma, dt), 0.0);
                 break;

        case 1 : payoff = std::max(strike - diff_simu_min_price(steps, stock, inRate, sigma, dt), 0.0);
    }

    return payoff;
}

double one_step (const double &old, double &inRate, double &vol, double &dt){
// the structure of the code, creat obj each time, is this too messy/expensive?
    std::normal_distribution<double> norDist(0.0, 1.0);
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    double one_draw = norDist(generator);
    double new_price = old * (1 + inRate * dt + vol * sqrt(dt) * one_draw);
    return new_price;
}

double diff_simu_min_price(const int &steps, double &start, double &inRate,
        double &vol, double &dt) {
    double minPrice = 100;
    double cum_lnS = 0;
    double current_price = start;
    std::normal_distribution<double> norDist(0.0, 1.0);
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    for (int count = 1; count <= steps;  count++) {
        double one_draw = norDist(generator);
        double dln_S = (inRate - 0.5 * std::pow(vol, 2)) * dt + vol * sqrt(dt) * one_draw;
        cum_lnS += dln_S;
        current_price = start * exp(cum_lnS);
        if (current_price < minPrice) {
            minPrice = current_price;
        }
    }

    return minPrice;
}

double end_price (const int &steps, double &start, double &inRate, double &vol, double &dt) {
    double current_price = start;
    for (int n = 1; n <= steps; n++) {
        current_price = one_step(current_price, inRate, vol, dt);
    }
    return current_price;
}



double max_price (const int &steps, double &start, double &inRate, double &vol, double &dt) {
    double maxPrice = -100;
    double current_price = start;
    for (int n = 1; n <= steps; n++) {
        current_price = one_step(current_price, inRate, vol, dt);
        if(current_price > maxPrice) {
            maxPrice = current_price;
        }
    }

    return maxPrice;
}

double min_price (const int &steps, double &start, double &inRate, double &vol, double &dt) {
    double minPrice = 1000;
    double current_price = start;
    for (int n = 1; n <= steps; n++) {
        current_price = one_step(current_price, inRate, vol, dt);
        if(current_price < minPrice) {
            minPrice = current_price;
        }
    }

    return minPrice;
}

