#include "cell.hpp"

#include <iostream>

int main() {
    std::cout << "sizeof(Cell): " << sizeof(tumopp::Cell) << "\n";
    tumopp::Cell cell({{1, 2, 3}}, 42);
    std::cout << tumopp::Cell::header() << "\n"
              << cell << "\n";
    return 0;
}
