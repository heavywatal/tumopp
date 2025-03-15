#include "cell.hpp"

#include <iostream>
#include <type_traits>

static_assert(std::is_nothrow_copy_constructible_v<tumopp::Cell>, "");
static_assert(std::is_nothrow_move_constructible_v<tumopp::Cell>, "");

int main() {
    std::cout << "sizeof(Cell): " << sizeof(tumopp::Cell) << "\n";
    tumopp::Cell cell({{1, 2, 3}}, 42);
    std::cout << tumopp::Cell::header() << "\n"
              << cell << "\n";
    return 0;
}
