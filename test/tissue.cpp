#include "tissue.hpp"

#include <iostream>

int main() {
    std::cout.precision(15);

    tumopp::Tissue tissue;
    tissue.grow(10);
    std::cout << tissue << "\n";
    tissue.write_history(std::cout);
    return 0;
}
