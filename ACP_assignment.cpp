#include <chrono>
#include <cmath>
#include <iostream>
#include <random>
#include <vector>

struct Option {
    double inRate, vol, expire, stock, strike;
    int steps;
    double dt, discount_factor;
};

double get_payoff(const int &gap, const bool &floating, const bool &call,
            const bool &arithmetic, const bool &continous, const Option &option);
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
    double running_avg = 0;
    double payoff = -1;
    std::vector<bool> argBool {true, false};
    std::vector<int> gaps {5, 22};
    for (int gap : gaps) {
        for (bool flo : argBool) {
            for (bool op_type : argBool) {
                for (bool ave_type : argBool) {
                    running_avg = 0;
                    for (int n = 1; n <= npaths; n++) {
                        payoff = get_payoff(gap, flo, op_type, ave_type, false, option);
                        running_avg = (running_avg * (n - 1) + payoff) / n;
                    }

                    double price = option.discount_factor * running_avg;
                    std::cout << "Option info : ";
                    std::cout << "Discrete sampling (each " << gap << " days), ";
                    std::cout << (flo ? "Floating strike, " : "Fixed strike, ");
                    std::cout << (ave_type ? "Arithmetic Average, " : "Geometric Average, ");
                    std::cout << (op_type ? "Call option" : "Put option") << std::endl;
                    std::cout << "Price : " << price << std::endl;
                }
            }
        }
    }

    for (bool flo : argBool) {
        for (bool op_type : argBool) {
            for (bool ave_type : argBool) {
                running_avg = 0;
                for (int n = 1; n <= npaths; n++) {
                    payoff = get_payoff(5, flo, op_type, ave_type, true, option);
                    running_avg = (running_avg * (n - 1) + payoff) / n;
                }

                double price = option.discount_factor * running_avg;
                std::cout << "Option info : ";
                std::cout << "Continuous sampling, ";
                std::cout << (flo ? "Floating strike, " : "Fixed strike, ");
                std::cout << (ave_type ? "Arithmetic Average, " : "Geometric Average, ");
                std::cout << (op_type ? "Call option" : "Put option") << std::endl;
                std::cout << "Price : " << price << std::endl;
            }
        }
    }

    return 0;
}

double get_payoff(const int &gap, const bool &floating, const bool &call, const bool &arithmetic, const bool &continous, const Option &option) {
    double payoff = -1;
    double old_price = option.stock;
    double new_price = one_step(old_price, option);
    double average = new_price;
    double tmp;
    old_price = new_price;
    for (int i = 2; i <= option.steps; i++) {
        new_price = one_step(old_price, option);
        if (continous){
            if (arithmetic){
                average = (average * (i - 1) + new_price) / i;
            } else {
                average = std::pow(average, (i - 1.0) / i) * std::pow(new_price, 1.0/i);
            }
        } else {
            int count = 0;
            if ((i - 1) % gap == 0) {
                count = (i - 1) / gap;
                if (arithmetic){
                    average = (average * count + new_price) / (count + 1);
                } else {
                    average = std::pow((std::pow(average, count) * new_price), 1.0/(count + 1));
                }
            }
        }

        old_price = new_price;
    }

    if (floating) {
        tmp = call ? new_price - average : average - new_price;
    } else {
        tmp = call ? average - option.strike : option.strike - average;
    }
    payoff = std::max(tmp, 0.0);

    return payoff;
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

