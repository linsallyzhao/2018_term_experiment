#include <iostream>
#include <cstdlib>
#include <ctime>

int main() {
    unsigned seed;
    std::cout << RAND_MAX << " is rand_max.\n";

    int total_n = 10000;
    int count = 0;
    seed = time(0);
    srand(seed);
    for (int n = 1; n <= total_n; n++){
        int num = rand() % 100;
        if (num < 3) {
            count++;
        }

    }
    std::cout << "Total count at last is " << count << std::endl;
    double ratio = static_cast<double>(count) / total_n;
    std::cout << ratio << " out of " << total_n << " was head. \n";

    return 0;
}


