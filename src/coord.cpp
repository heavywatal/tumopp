#include "coord.hpp"

#include <cmath>
#include <numeric>
#include <algorithm>
#include <bitset>
#include <stdexcept>

namespace tumopp {

namespace {

//! constexpr pow function for integers
constexpr size_t ipow(size_t base, size_t exponent) noexcept {
  return exponent <= 0u ? 1u
       : exponent == 1u ? base
       : base * ipow(base, --exponent);
}

template <class T> inline
T sum(const std::array<T, MAX_DIM>& v) {
    return std::accumulate(v.begin(), v.end(), T{});
}

template <class T> inline
T max(const std::array<T, MAX_DIM>& v) {
    return *std::max_element(v.begin(), v.end());
}

template <class T> inline
std::array<T, MAX_DIM> abs(const std::array<T, MAX_DIM>& v) {
    std::array<T, MAX_DIM> out{};
    std::transform(v.begin(), v.end(), out.begin(),
                   [](T x){return std::abs(x);});
    return out;
}

//! Euclidean distance
template <class T> inline
double _euclidean_distance(const std::array<T, MAX_DIM>& v) {
    return std::sqrt(sum((v * v)));
}

}// namespace

Coord::Coord(unsigned d): dimensions_(d) {
    if (d < 2U || 3U < d) {
        throw std::runtime_error("Invalid value for dimensions");
    }
}

Neumann::Neumann(const unsigned d): Coord(d) {
    directions_.reserve(2U * dimensions_);
    for (unsigned i = 0; i < dimensions_; ++i) {
        coord_t v{};
        v[i] = 1;
        directions_.push_back(v);
        v[i] = -1;
        directions_.push_back(v);
    }
    auto dirmax = static_cast<unsigned>(directions_.size()) - 1;
    dist_direction_.param(decltype(dist_direction_)::param_type(0, dirmax));
}

Moore::Moore(const unsigned d): Coord(d) {
    directions_.reserve(ipow(3u, dimensions_) - 1u);
    for (const int x: {-1, 0, 1}) {
        for (const int y: {-1, 0, 1}) {
            if (dimensions_ == 2U) {
                if (x == 0 && y == 0) continue;
                directions_.push_back({{x, y, 0}});
                continue;
            }
            for (const int z: {-1, 0, 1}) {
                if (x == 0 && y == 0 && z == 0) continue;
                directions_.push_back({{x, y, z}});
            }
        }
    }
    auto dirmax = static_cast<unsigned>(directions_.size()) - 1;
    dist_direction_.param(decltype(dist_direction_)::param_type(0, dirmax));
}

Hexagonal::Hexagonal(const unsigned d): Coord(d) {
    coord_t v{{-1, 0, 1}};
    directions_.reserve(6 * (dimensions_ - 1));
    do {
        directions_.push_back({{v[0], v[1], 0}});
    } while (std::next_permutation(std::begin(v), std::end(v)));
    if (dimensions_ > 2u) {
        directions_.push_back({{0, 0, -1}});
        directions_.push_back({{1, 0, -1}});
        directions_.push_back({{1, -1, -1}});
        directions_.push_back({{0, 0, 1}});
        directions_.push_back({{-1, 0, 1}});
        directions_.push_back({{-1, 1, 1}});
    }
    auto dirmax = static_cast<unsigned>(directions_.size()) - 1;
    dist_direction_.param(decltype(dist_direction_)::param_type(0, dirmax));
}

int Neumann::graph_distance(const coord_t& v) const {
    return sum(abs(v));
}

int Moore::graph_distance(const coord_t& v) const {
    return max(abs(v));
}

std::array<double, MAX_DIM> Coord::continuous(const coord_t& v) const {
    std::array<double, MAX_DIM> out{};
    std::copy(std::begin(v), std::end(v), std::begin(out));
    return out;
}

std::array<double, MAX_DIM> Hexagonal::continuous(const coord_t& v) const {
    std::array<double, MAX_DIM> true_pos{};
    true_pos[0] += static_cast<double>(v[0]);
    true_pos[1] += static_cast<double>(v[1]);
    true_pos[1] += true_pos[0] * 0.5;
    true_pos[0] *= std::sqrt(3.0 / 4.0);
    if (dimensions_ > 2U) {
        true_pos[2] += static_cast<double>(v[2]);
        true_pos[0] += true_pos[2] / std::sqrt(3.0);
        true_pos[2] *= std::sqrt(2.0 / 3.0);
    }
    return true_pos;
}

int Hexagonal::graph_distance(const coord_t& v) const {
    int d = std::max(max(abs(v)), std::abs(v[0] + v[1]));
    if (dimensions_ > 2U) {
        return std::max(d, std::abs(v[0] + v[2]));
    }
    return d;
}

double Coord::euclidean_distance(const coord_t& v) const {
    return _euclidean_distance(v);
}

double Hexagonal::euclidean_distance(const coord_t& v) const {
    return _euclidean_distance(continuous(v));
}

std::vector<coord_t> Coord::core() const {
    const size_t n = ipow(2u, dimensions_);
    std::vector<coord_t> output;
    output.reserve(n);
    for (size_t i=0; i<n; ++i) {
        std::bitset<MAX_DIM> bs(i);
        coord_t v{};
        for (size_t j=0; j<dimensions_; ++j) {
            v[j] = static_cast<int>(bs[j]);
        }
        output.push_back(v);
    }
    return output;
}

std::vector<coord_t> Hexagonal::core() const {
    std::vector<coord_t> output = Neumann(dimensions_).core();
    if (dimensions_ == 3U) {
        output.resize(3);
        output.push_back({{1, 0, -1}});
    }
    return output;
}

std::vector<coord_t> Coord::sphere(const size_t n) const {
    std::vector<coord_t> output;
    if (dimensions_ == 2U) {
        const int lim = 9;
        // radius 9: regular: 253, hex: 281
        output.reserve(281);
        for (int x=-lim; x<=lim; ++x) {
            for (int y=-lim; y<=lim; ++y) {
                coord_t v{{x, y, 0}};
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
                    coord_t v{{x, y, z}};
                    if (euclidean_distance(v) <= lim) {
                        output.push_back(v);
                    }
                }
            }
        }
    }
    std::sort(output.begin(), output.end(),
        [this](const coord_t& lhs, const coord_t& rhs){
            return euclidean_distance(lhs) < euclidean_distance(rhs);
    });
    output.resize(n);
    return output;
}

coord_t Coord::outward(const coord_t& v) const {
    const auto candidates = neighbors(v);
    return *std::max_element(candidates.begin(), candidates.end(),
                             [this](const coord_t& lhs, const coord_t& rhs) {
        return euclidean_distance(lhs) < euclidean_distance(rhs);
    });
}

} // namespace tumopp
