#include "tissue.hpp"

#include <wtl/iostr.hpp>

//! Test coordinate system
template <class T> inline
void test_coordinate(const std::valarray<int>& v) {
    T coord(static_cast<unsigned>(v.size()));
    std::cout << coord.core() << "\n";
    std::cout << coord.sphere(20) << "\n";
    for (auto x: coord.directions()) {
        std::cout << x  << ": " << coord.euclidean_distance(x) << "\n";
    }
    std::cout << coord.neighbors(v) << "\n";
    std::cout << coord.graph_distance(v) << "\n";
    std::cout << coord.euclidean_distance(v) << "\n";
    std::cout << coord.radius(100) << "\n";
    std::cout << coord.radius(1000) << "\n";
}

int main() {
    std::cout.precision(15);
    tumopp::Cell::init_distributions();

    tumopp::Tissue tissue;
    tissue.grow(10);
    std::cout << tissue << "\n";
    std::cout << tissue.specimens();

    const std::valarray<int> v2{3, -2};
    test_coordinate<tumopp::Neumann>(v2);
    test_coordinate<tumopp::Moore>(v2);
    test_coordinate<tumopp::Hexagonal>(v2);

    const std::valarray<int> v3{3, -2, 1};
    test_coordinate<tumopp::Neumann>(v3);
    test_coordinate<tumopp::Moore>(v3);
    test_coordinate<tumopp::Hexagonal>(v3);
    return 0;
}
