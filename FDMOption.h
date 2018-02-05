#include <iostream>
#include <vector>
#include <cmath>

#ifndef FDMOPTION_H
#define FDMOPTION_H

class FDMOption{
    private:
        double beginningStockPrice, strike, term;
        bool call;
        vector<double> a_n, b_n, c_n, A_n, B_n, C_n;

    public:
        FDMOption(double S0, double E, double T, bool C) {
            beginningStockPrice = S0;
            strike = E;
            term = T;
            call = C;
        }

        //No default number for any private parameter since there is no
        //reasonable default value for any. When you create a option you have
        //to define these otherwise it is not a option

        //For the same reason,  no function to change these four after a option
        //is created since you change any of these you change the option

        double getBeginningStockPrice() const{
            return beginningStockPrice;
        }

        double getStrikePrice() const{
            return strike;
        }

        double getExpiry() const{
            return term;
        }

        void getOptionType() const{
            std::cout << "This is a European ";
            if (call)
                std::cout << "call ";
            else
                std::cout << "put ";
            std::cout << "option." << std::endl;
        }

        void fill_vectors (double riskFree, double vol, double div, double weight, int numTimeStep) {
            double deltaT = term / numTimeStep;
            //Don't like the vector method. It is my first thought because it
            //is easy to handle but it is not efficient. I want to dynamically
            //allocate six arraies here and I will know how long they need to
            //be since I know numTimeStep. And at the end, to pass this info to
            //the pricing function, I can pass a array contains six
            //addresses. Do I also need to pass how long each one is?
            //I do think it is a good way to make these arraies member variable
            //of class so any other function can use them with being passed.
            //But can a mumber function creat new member variables?
            b_n.push_back(1 + weight * riskFree * deltaT);
            c_n.push_back(0);
            B_n.push_back(1 - (1 - weight) * riskFree * deltaT);
            C_n.push_back(0);
            //There should be some safety guard for a_n and A_n, you should never
            //use a_n[0] or A_n[0], neither c_n[numTimeStep] or
            //C_n[numTimeSteps], but for C_n and c_n it won't be that long

            for (int index = 1; index < numTimeStep; index++){
                double block_1 = (pow(vol, 2) * pow(index, 2);
                double block_2 = index * ( riskFree - div);
                a_n.push_back(-0.5 * weight * (block_1 - block_2) * deltaT);
                b_n.push_back(1 + weight * (block_1 + riskFree) * deltaT);
                c_n.push_back(-0.5 * weight * (block_1 + block_2) * deltaT);
                A_n.push_back(0.5 * (1 - weight) * (block_1 - block_2) * deltaT);
                B_n.push_back(1 - (1 - weight) * (block_1 + riskFree) * deltaT);
                C_n.push_back(0.5 * (1 - weight) * (block_1 + block_2) * deltaT);
            }
            //Need boundary conditon here
            //This whole function has not been tested yet because I don't know
            //how to return things

        }



};

#endif


