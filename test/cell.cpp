#include "cell.hpp"

#include <iostream>

int main() {
    tumopp::Cell cell({1, 2, 3});
    std::cout << tumopp::Cell::header() << "\n"
              << cell << std::endl;
}
