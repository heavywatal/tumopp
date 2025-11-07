#include "tissue.hpp"

#include <fmt/format.h>

#include <iostream>

int main() {
    tumopp::Tissue tissue;
    tissue.grow(10);
    fmt::println("{}", tissue);
    tissue.write_history(std::cout);
    return 0;
}
