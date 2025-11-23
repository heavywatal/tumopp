#include "coord.hpp"

#include <fmt/format.h>
#include <fmt/ranges.h>
#include <typeinfo>

template <class T> inline
void test_coordinate(int dim, const tumopp::coord_t& v) {
    fmt::println("{} {}D ################", typeid(T).name(), dim);
    T coord_func(dim);
    for (auto c: coord_func.directions()) {
        fmt::println("{}: {:.9g}", c, coord_func.euclidean_distance(c));
    }
    fmt::println("core:        {}", coord_func.core());
    fmt::println("dist_g:      {}", coord_func.graph_distance(v));
    fmt::println("dist_e:      {:.9g}", coord_func.euclidean_distance(v));
}

void test_dimension(int dim, tumopp::coord_t&& v) {
    test_coordinate<tumopp::Neumann>(dim, v);
    test_coordinate<tumopp::Moore>(dim, v);
    test_coordinate<tumopp::Hexagonal>(dim, v);
}

int main() {
    test_dimension(1, {1, 0, 0});
    test_dimension(2, {1, -2, 0});
    test_dimension(3, {1, -2, 3});
    return 0;
}
