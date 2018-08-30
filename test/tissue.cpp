#include "tissue.hpp"

#include <iostream>

int main() {
    std::cout.precision(15);

    tumopp::Tissue tissue;
    tissue.grow(10);
    std::cout << tissue << "\n";
    std::cout << tissue.history().rdbuf();
    return 0;
}
