#include "tissue.hpp"

#include <wtl/iostr.hpp>

#include <iostream>

//! Test coordinate system
template <class T> inline
void test_coordinate(const std::valarray<int>& v) {
    T coord(static_cast<unsigned int>(v.size()));
    std::cout << coord.core() << std::endl;
    std::cout << coord.sphere(20) << std::endl;
    for (auto x: coord.directions()) {
        std::cout << x  << ": " << coord.euclidean_distance(x) << std::endl;
    }
    std::cout << coord.neighbors(v) << std::endl;
    std::cout << coord.graph_distance(v) << std::endl;
    std::cout << coord.euclidean_distance(v) << std::endl;
    std::cout << coord.radius(100) << std::endl;
    std::cout << coord.radius(1000) << std::endl;
}

int main() {
    std::cout.precision(15);
    tumopp::Cell::init_distributions();

    tumopp::Tissue tissue;
    tissue.grow(10);
    std::cout << tissue << std::endl;
    std::cout << tissue.header();
    tissue.write(std::cout);

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
