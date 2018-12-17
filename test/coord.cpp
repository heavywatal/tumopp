#include "coord.hpp"

#include <wtl/iostr.hpp>
#include <typeinfo>

template <class T> inline
void test_coordinate(unsigned dim, const tumopp::coord_t& v) {
    std::cout << typeid(T).name() << " " << dim << "D ################\n";
    T coord_func(dim);
    for (auto x: coord_func.directions()) {
        std::cout << x << ": " << coord_func.euclidean_distance(x) << "\n";
    }
    std::cout << "core:        " << coord_func.core() << "\n";
    std::cout << "neighbors:   " << coord_func.neighbors(v) << "\n";
    std::cout << "dist_g:      " << coord_func.graph_distance(v) << "\n";
    std::cout << "dist_e:      " << coord_func.euclidean_distance(v) << "\n";
}

void test_dimension(unsigned dim, tumopp::coord_t&& v) {
    test_coordinate<tumopp::Neumann>(dim, v);
    test_coordinate<tumopp::Moore>(dim, v);
    test_coordinate<tumopp::Hexagonal>(dim, v);
}

int main() {
    std::cout.precision(9);
    test_dimension(1u, {1, 0, 0});
    test_dimension(2u, {1, -2, 0});
    test_dimension(3u, {1, -2, 3});
    return 0;
}
