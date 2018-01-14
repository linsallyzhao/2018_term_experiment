#include <chrono>
#include <cmath>
#include <iostream>
#include <random>

struct Option {
    double inRate, vol, expire, stock, strike;
    int steps;
    double dt, discount_factor;
};

double get_payoff(const int &type, const Option &put);
double one_step(const double &old, const Option &put);
double diff_simu_min_price(const Option &put);
double end_price(const Option &put);
double min_price(const Option &put);
double max_price(const Option &put);

int main() {
    // swear I am going to put these into a strucure
    Option put;
    put.inRate = 0.1;
    put.vol = 0.3;
    put.expire = 0.5;
    put.stock = 1;
    put.strike = 1;
    put.steps = 50;
    put.dt = put.expire / put.steps;
    put.discount_factor = exp(-put.inRate * put.expire);
    int npaths = 1000000;
    double running_avg = 0;
    double payoff;
    for (int type = 0; type <= 1; type++) {
        for (int n = 1; n <= npaths; n++) {
            payoff = get_payoff(type, put);
            running_avg = (running_avg * (n - 1) + payoff) / n;
        }

        double price = put.discount_factor * running_avg;
        if (type == 0) {
            std::cout << "Monte carlo price for lookback min_price put is "
                      << price << std::endl;
        } else if (type == 1) {
            std::cout << "Monte carlo price for lookback diff_simu put is "
                      << price << std::endl;
        }
    }

    return 0;
}

double get_payoff(const int &type, const Option &put) {
    double payoff;

    switch (type) {
        case 0:
            payoff = std::max(put.strike - min_price(put), 0.0);
            break;

        case 1:
            payoff = std::max(put.strike - diff_simu_min_price(put), 0.0);
            break;
    }

    return payoff;
}

double one_step(const double &old, const Option &put) {
    std::normal_distribution<double> norDist(0.0, 1.0);
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    double one_draw = norDist(generator);
    double new_price =
        old * (1 + put.inRate * put.dt + put.vol * sqrt(put.dt) * one_draw);
    return new_price;
}

double diff_simu_min_price(const Option &put) {
    double minPrice = 100;
    double cum_lnS = 0;
    double current_price = put.stock;
    std::normal_distribution<double> norDist(0.0, 1.0);
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    for (int count = 1; count <= put.steps; count++) {
        double one_draw = norDist(generator);
        double dln_S = (put.inRate - 0.5 * std::pow(put.vol, 2)) * put.dt +
                       put.vol * sqrt(put.dt) * one_draw;
        cum_lnS += dln_S;
        current_price = put.stock * exp(cum_lnS);
        if (current_price < minPrice) {
            minPrice = current_price;
        }
    }

    return minPrice;
}

double end_price(const Option &put) {
    double current_price = put.stock;
    for (int n = 1; n <= put.steps; n++) {
        current_price = one_step(current_price, put);
    }
    return current_price;
}

double max_price(const Option &put) {
    double maxPrice = -100;
    double current_price = put.stock;
    for (int n = 1; n <= put.steps; n++) {
        current_price = one_step(current_price, put);
        if (current_price > maxPrice) {
            maxPrice = current_price;
        }
    }

    return maxPrice;
}

double min_price(const Option &put) {
    double minPrice = 1000;
    double current_price = put.stock;
    for (int n = 1; n <= put.steps; n++) {
        current_price = one_step(current_price, put);
        if (current_price < minPrice) {
            minPrice = current_price;
        }
    }

    return minPrice;
}
