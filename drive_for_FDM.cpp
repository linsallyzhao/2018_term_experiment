#include "FDMOption.h"

int main() {
    FDMOption option(100, 105, 1.0, 1);
    
    std::cout << "The beginning stock price is " << option.getBeginningStockPrice() << std::endl;
    std::cout << "The strike price is " << option.getStrikePrice() << std::endl;
    std::cout << "The term of the option is " << option.getExpiry() << std::endl;

    option.getOptionType();


}
