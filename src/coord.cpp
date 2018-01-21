#include "coord.hpp"

#include <boost/math/constants/constants.hpp>

#include <cmath>
#include <cstdlib> // std::abs for int
#include <iostream>
#include <bitset>

namespace tumopp {

//! Mathematical constant
constexpr double PI = boost::math::constants::pi<double>();

//! constexpr pow function for integers
constexpr size_t ipow(size_t base, size_t exponent) noexcept {
  return exponent <= 0u ? 1u
       : exponent == 1u ? base
       : base * ipow(base, --exponent);
}

Neumann::Neumann(const unsigned int d): Coord(d) {
    directions_.reserve(2U * dimensions_);
    std::valarray<int> v(dimensions_);
    v[v.size() - 1] += 1;
    do {
        directions_.push_back(v);
    } while (std::next_permutation(std::begin(v), std::end(v)));
    v = 0;
    v[0] -= 1;
    do {
        directions_.push_back(v);
    } while (std::next_permutation(std::begin(v), std::end(v)));
    auto dirmax = static_cast<unsigned int>(directions_.size()) - 1;
    dist_direction_.param(decltype(dist_direction_)::param_type(0, dirmax));
}

Moore::Moore(const unsigned int d): Coord(d) {
    directions_.reserve(ipow(3u, dimensions_) - 1u);
    for (const int x: {-1, 0, 1}) {
        for (const int y: {-1, 0, 1}) {
            if (dimensions_ == 2U) {
                if (x == 0 && y == 0) continue;
                directions_.push_back({x, y});
                continue;
            }
            for (const int z: {-1, 0, 1}) {
                if (x == 0 && y == 0 && z == 0) continue;
                directions_.push_back({x, y, z});
            }
        }
    }
    auto dirmax = static_cast<unsigned int>(directions_.size()) - 1;
    dist_direction_.param(decltype(dist_direction_)::param_type(0, dirmax));
}

Hexagonal::Hexagonal(const unsigned int d): Coord(d) {
    std::valarray<int> v{-1, 0, 1};
    directions_.reserve(6 * (dimensions_ - 1));
    if (dimensions_ == 2U) {
        do {
            directions_.push_back({v[0], v[1]});
        } while (std::next_permutation(std::begin(v), std::end(v)));
    }
    else {
        do {
            directions_.push_back({v[0], v[1], 0});
        } while (std::next_permutation(std::begin(v), std::end(v)));
        directions_.push_back({0, 0, -1});
        directions_.push_back({1, 0, -1});
        directions_.push_back({1, -1, -1});
        directions_.push_back({0, 0, 1});
        directions_.push_back({-1, 0, 1});
        directions_.push_back({-1, 1, 1});
    }
    auto dirmax = static_cast<unsigned int>(directions_.size()) - 1;
    dist_direction_.param(decltype(dist_direction_)::param_type(0, dirmax));
}

int Neumann::graph_distance(const std::valarray<int>& v) const {
    const std::valarray<int> absv = std::abs(v);
    return absv.sum();
}

int Moore::graph_distance(const std::valarray<int>& v) const {
    return std::abs(v).max();
}

int Hexagonal::graph_distance(const std::valarray<int>& v) const {
    int d = std::max(std::abs(v).max(), std::abs(v[0] + v[1]));
    if (dimensions_ > 2U) {
        return std::max(d, std::abs(v[0] + v[2]));
    }
    return d;
}

double Hexagonal::euclidean_distance(const std::valarray<int>& v) const {
    std::valarray<double> true_pos(dimensions_);
    true_pos[0] += static_cast<double>(v[0]);
    true_pos[1] += static_cast<double>(v[1]);
    true_pos[1] += true_pos[0] * 0.5;
    true_pos[0] *= std::sqrt(3.0 / 4.0);
    if (dimensions_ > 2U) {
        true_pos[2] += static_cast<double>(v[2]);
        true_pos[0] += true_pos[2] / sqrt(3.0);
        true_pos[2] *= std::sqrt(2.0 / 3.0);
    }
    return _euclidean_distance(true_pos);
}

std::vector<std::valarray<int>> Coord::core() const {
    const size_t n = ipow(2u, dimensions_);
    std::vector<std::valarray<int>> output;
    output.reserve(n);
    for (size_t i=0; i<n; ++i) {
        std::bitset<3> bs(i);
        std::valarray<int> v(dimensions_);
        for (size_t j=0; j<dimensions_; ++j) {
            v[j] = static_cast<int>(bs[j]);
        }
        output.push_back(v);
    }
    return output;
}

std::vector<std::valarray<int>> Hexagonal::core() const {
    std::vector<std::valarray<int>> output = Neumann(dimensions_).core();
    if (dimensions_ == 3U) {
        output.resize(3);
        output.push_back({1, 0, -1});
    }
    return output;
}

double Coord::radius(const size_t nodes) const {
    double x = nodes;
    x /= PI;
    if (dimensions_ == 2U) {
        // S = pi r^2
        return std::sqrt(x);
    } else {
        // V = 4/3 pi r^3
        return std::pow(x *= (3.0 / 4.0), 1.0 / 3.0);
    }
}

double Hexagonal::radius(const size_t nodes) const {
    const static double rate2d = std::pow(std::sqrt(3.0 / 4.0), 1.0 / 3.0);
    const static double rate3d = std::pow(std::sqrt(0.5), 1.0 / 3.0);
    if (dimensions_ == 2U) {
        return Coord::radius(nodes) * rate2d;
    } else {
        return Coord::radius(nodes) * rate3d;
    }
}

double Coord::cross_section(size_t nodes) const {
    double r = radius(nodes);
    return PI * (r *= r);
}

std::vector<std::valarray<int>> Coord::sphere(const size_t n) const {
    std::vector<std::valarray<int>> output;
    if (dimensions_ == 2U) {
        const int lim = 9;
        // radius 9: regular: 253, hex: 281
        output.reserve(281);
        for (int x=-lim; x<=lim; ++x) {
            for (int y=-lim; y<=lim; ++y) {
                std::valarray<int> v = {x, y};
                if (euclidean_distance(v) <= lim) {
                    output.push_back(v);
                }
            }
        }
    } else {
        const int lim = 4;
        // radius 4: regular: 257, hex: 357
        output.reserve(357);
        for (int x=-lim; x<=lim; ++x) {
            for (int y=-lim; y<=lim; ++y) {
                for (int z=-lim; z<=lim; ++z) {
                    std::valarray<int> v = {x, y, z};
                    if (euclidean_distance(v) <= lim) {
                        output.push_back(v);
                    }
                }
            }
        }
    }
    std::sort(output.begin(), output.end(),
        [this](const std::valarray<int>& lhs, const std::valarray<int>& rhs){
            return euclidean_distance(lhs) < euclidean_distance(rhs);
    });
    output.resize(n);
    return output;
}

}
