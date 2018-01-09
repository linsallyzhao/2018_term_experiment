#include <iostream>
#include <random>
#include <cmath>
#include <vector>
#include <numeric>
#include <algorithm>
#include <chrono>
using namespace std;

// want to change these to const function. They don't modify things
double one_step (double &old, double &inRate, double &vol, double &dt);
double end_price (int &steps, double &start, double &inRate, double &vol, double &dt);
//double max_price ();


int main() {
    double inRate = 0.05, sigma = 0.2, expire = 1;
    double stock = 100, strike = 100;
    int steps = 252; 
    double dt = expire / steps;
    double discount_factor = exp(-inRate * expire);
    int npaths = 10000;
    vector<double> pay_offs;
    double payoff;
    for (int n = 1; n <= npaths; n++) {
        payoff = max(end_price(steps, stock, inRate, sigma, dt) - strike, 0.0);
//cout << "pay off " << payoff << endl;
        pay_offs.push_back(payoff);
    }

    double avg = accumulate(pay_offs.begin(), pay_offs.end(), 0.0) / pay_offs.size();
    double price = discount_factor * avg;

    cout << "Monte carlo price is " << price << endl;

    return 0;
}


double one_step (double &old, double &inRate, double &vol, double &dt){
// the structure of the code, creat obj each time, is this too messy/expensive?
    //minstd_rand0 generator;
    normal_distribution<double> norDist(0.0, 1.0);
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator(seed);
    double one_draw = norDist(generator);
//cout << "One draw " << one_draw << endl;
    double new_price = old * (1 + inRate * dt + vol * sqrt(dt) * one_draw);
//cout << "update price, new price is " << new_price << endl;
    return new_price;
}

double end_price (int &steps, double &start, double &inRate, double &vol, double &dt) {
    double current_price = start;
    for (int n = 1; n <= steps; n++) {
        current_price = one_step(current_price, inRate, vol, dt);
    }
//cout << "end price is " << current_price << endl;
    return current_price;
}



