#include "tissue.hpp"

#include <fmt/format.h>

#include <iostream>

int main() {
    tumopp::Tissue tissue;
    tissue.grow(10);
    fmt::println("{}", tissue);
    std::cout << tissue.str_history();
    return 0;
}
