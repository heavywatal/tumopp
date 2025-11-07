#include "cell.hpp"

#include <fmt/format.h>

#include <type_traits>

static_assert(std::is_nothrow_copy_constructible_v<tumopp::Cell>, "");
static_assert(std::is_nothrow_move_constructible_v<tumopp::Cell>, "");

int main() {
    fmt::println("sizeof(Cell): {}", sizeof(tumopp::Cell));
    tumopp::Cell cell({{1, 2, 3}}, 42);
    fmt::print("{}\n{}\n", tumopp::Cell::header(), cell);
    return 0;
}
