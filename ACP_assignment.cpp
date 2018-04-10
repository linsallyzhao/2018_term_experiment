#include <chrono>
#include <cmath>
#include <iostream>
#include <random>
#include <vector>
#include <unordered_map>
#include <string>
struct Option {
    double inRate, vol, expire, stock, strike;
    int steps;
    double dt, discount_factor;
};

std::unordered_map<std::string, double> get_payoffs(const int &gap, const Option &option);
double one_step(const double &old, const Option &option);

int main() {
    Option option;
    option.inRate = 0.05;
    option.vol = 0.2;
    option.expire = 1;
    option.stock = 100;
    option.strike = 100;
    option.steps = 252;
    option.dt = option.expire / option.steps;
    option.discount_factor = exp(-option.inRate * option.expire);
    int npaths = 100000;
    std::unordered_map<std::string, double> payoffs;
    std::vector<int> gaps {5, 22};
    std::unordered_map<std::string, double> running_avgs;
    std::unordered_map<std::string, double> running_sqr_avgs;
    int gap = 5;

    for (const auto & [ key, value ] : payoffs) {
        running_avgs[key] = 0.0;
        running_sqr_avgs[key]  = 0.0;
    }

    for (int n = 1; n <= npaths; n++) {
        payoffs = get_payoffs(gap, option);
        for (const auto & [ key, value ] : payoffs) {
            running_avgs[key] = (running_avgs[key] * (n - 1) + value) / n;
            running_sqr_avgs[key] = (running_sqr_avgs[key] * (n - 1) + value * value) / n;
        }
    }

    double price = 0.0;
    for (const auto & [ key, value ] : running_avgs) {
        price = option.discount_factor * value;
        std::cout << key << " " << price << std::endl;
        double std_payoffs = std::pow(running_sqr_avgs[key] *
                    std::pow(option.discount_factor, 2) - std::pow(price, 2), 0.5);
        double std_error = std_payoffs / std::pow(npaths, 0.5);

        std::cout << "standard error is " << std_error << std::endl;

    }


    return 0;
}


std::unordered_map<std::string, double> get_payoffs(const int &gap, const Option &option) {
    double payoff = -1;
    double old_price = option.stock;
    double new_price = one_step(old_price, option);
    std::unordered_map<std::string, double> averages = {
        {"con_ari", new_price},
        {"con_geo", new_price},
        {"dis_ari", new_price},
        {"dis_geo", new_price}
    };
    double tmp;
    old_price = new_price;
    for (int i = 2; i <= option.steps; i++) {
        new_price = one_step(old_price, option);
        averages["con_ari"] = (averages["con_ari"] * (i - 1) + new_price) / i;
        averages["con_geo"] = std::pow(averages["con_geo"], (i - 1.0) / i) * std::pow(new_price, 1.0/i);

        int count = 0;
        if ((i - 1) % gap == 0) {
            count = (i - 1) / gap;
            averages["dis_ari"] = (averages["dis_ari"] * count + new_price) / (count + 1);
            averages["dis_geo"] = std::pow((std::pow(averages["dis_geo"], count) * new_price), 1.0/(count + 1));
        }
        old_price = new_price;
    }

    std::unordered_map<std::string, double> payoffs;
    payoffs["con_ari_flo_call"] = std::max(new_price - averages["con_ari"], 0.0);
    payoffs["con_ari_flo_put"] = std::max(averages["con_ari"] - new_price, 0.0);
    payoffs["con_ari_fix_call"] = std::max(averages["con_ari"] - option.strike, 0.0);
    payoffs["con_ari_fix_put"] = std::max(option.strike -  averages["con_ari"], 0.0);

    payoffs["con_geo_flo_call"] = std::max(new_price - averages["con_geo"], 0.0);
    payoffs["con_geo_flo_put"] = std::max(averages["con_geo"] - new_price, 0.0);
    payoffs["con_geo_fix_call"] = std::max(averages["con_geo"] - option.strike, 0.0);
    payoffs["con_geo_fix_put"] = std::max(option.strike -  averages["con_geo"], 0.0);

    payoffs["dis_ari_flo_call"] = std::max(new_price - averages["dis_ari"], 0.0);
    payoffs["dis_ari_flo_put"] = std::max(averages["dis_ari"] - new_price, 0.0);
    payoffs["dis_ari_fix_call"] = std::max(averages["dis_ari"] - option.strike, 0.0);
    payoffs["dis_ari_fix_put"] = std::max(option.strike -  averages["dis_ari"], 0.0);

    payoffs["dis_geo_flo_call"] = std::max(new_price - averages["dis_geo"], 0.0);
    payoffs["dis_geo_flo_put"] = std::max(averages["dis_geo"] - new_price, 0.0);
    payoffs["dis_geo_fix_call"] = std::max(averages["dis_geo"] - option.strike, 0.0);
    payoffs["dis_geo_fix_put"] = std::max(option.strike -  averages["dis_geo"], 0.0);
    return payoffs;
}

double one_step(const double &old, const Option &option) {
    std::normal_distribution<double> norDist(0.0, 1.0);
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    double one_draw = norDist(generator);
    double new_price =
        old * (1 + option.inRate * option.dt + option.vol * sqrt(option.dt) * one_draw);
    return new_price;
}

