#include "coord.hpp"

#include <wtl/iostr.hpp>

template <class T, unsigned N> inline
void test_coordinate(const tumopp::coord_t& v) {
    T coord_func(N);
    std::cout << "core:        " << coord_func.core() << "\n";
    std::cout << "sphere:      " << coord_func.sphere(20) << "\n";
    for (auto x: coord_func.directions()) {
        std::cout << x << ": " << coord_func.euclidean_distance(x) << "\n";
    }
    std::cout << "neighbors:   " << coord_func.neighbors(v) << "\n";
    std::cout << "dist_g:      " << coord_func.graph_distance(v) << "\n";
    std::cout << "dist_e:      " << coord_func.euclidean_distance(v) << "\n";
    std::cout << "radius(1e2): " << coord_func.radius(100) << "\n";
    std::cout << "radius(1e3): " << coord_func.radius(1000) << "\n";
}

int main() {
    std::cout.precision(9);
    const tumopp::coord_t v2{{3, -2, 0}};
    test_coordinate<tumopp::Neumann, 2>(v2);
    test_coordinate<tumopp::Moore, 2>(v2);
    test_coordinate<tumopp::Hexagonal, 2>(v2);

    const tumopp::coord_t v3{{3, -2, 1}};
    test_coordinate<tumopp::Neumann, 3>(v3);
    test_coordinate<tumopp::Moore, 3>(v3);
    test_coordinate<tumopp::Hexagonal, 3>(v3);
    return 0;
}
