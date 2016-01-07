// -*- mode: c++; coding: utf-8 -*-
/*! @file coord.hpp
    @brief Coordinate system
*/
#pragma once
#ifndef COORD_HPP_
#define COORD_HPP_

#include <cassert>
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>
#include <random>


/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

template <class T>
std::vector<T>& operator+=(std::vector<T>& lhs, const std::vector<T>& rhs) {
    assert(lhs.size() == rhs.size());
    std::transform(lhs.begin(), lhs.end(), rhs.begin(),
                   lhs.begin(), std::plus<T>());
    return lhs;
}

template <class T>
std::vector<T>& operator-=(std::vector<T>& lhs, const std::vector<T>& rhs) {
    assert(lhs.size() == rhs.size());
    std::transform(lhs.begin(), lhs.end(), rhs.begin(),
                   lhs.begin(), std::minus<T>());
    return lhs;
}

template <class T>
std::vector<T> operator+(std::vector<T> lhs, const std::vector<T>& rhs) {
    return lhs += rhs;
}

template <class T>
std::vector<T> operator-(std::vector<T> lhs, const std::vector<T>& rhs) {
    return lhs -= rhs;
}

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

class Coord {
  public:
    // methods
    const std::vector<std::vector<int>>& directions() const {return directions_;}
    size_t max_neighbors() const {return max_neighbors_;}

    std::vector<int> origin() const {
        return std::vector<int>(dimensions);
    }
    std::vector<std::vector<int>> neighbors(const std::vector<int>& v) const {
        std::vector<std::vector<int>> output = directions_;
        for (auto& d: output) {
            d += v;
        }
        return output;
    }
    template <class RNG> inline
    std::vector<int> random_neighbor(const std::vector<int>& v, RNG& rng) const {
        std::uniform_int_distribution<ptrdiff_t> uniform(0, max_neighbors_ - 1);
        return v + directions_[uniform(rng)];
    }
    std::vector<int> outward(const std::vector<int>& v) const {
        const auto candidates = neighbors(v);
        return *std::max_element(candidates.begin(), candidates.end(),
                                 [this](const std::vector<int>& lhs, const std::vector<int>& rhs) {
            return euclidean_distance(lhs) < euclidean_distance(rhs);
        });
    }

    // virtual methods
    virtual size_t graph_distance(const std::vector<int>& v) const = 0;
    virtual double euclidean_distance(const std::vector<int>& v) const {
        return _euclidean_distance(v);
    }
    virtual double radius(size_t volume) const {
        double x = volume;
        x /= M_PI;
        if (dimensions == 2) {
            // S = pi r^2
            return std::sqrt(x);
        } else {
            // V = 4/3 pi r^3
            x *= 3.0;
            x /= 4.0;
            return std::pow(x, 1.0/3.0);
        }
    }
    virtual std::vector<std::vector<int>> core() const {
        const size_t n = std::pow(2, dimensions);
        std::vector<std::vector<int>> output;
        output.reserve(n);
        for (size_t i=0; i<n; ++i) {
            std::bitset<3> bs(i);
            std::vector<int> v(dimensions);
            for (size_t j=0; j<dimensions; ++j) {
                v[j] = static_cast<int>(bs[j]);
            }
            output.push_back(v);
        }
        return output;
    }

    virtual ~Coord() = default;

    // properties
    const size_t dimensions;

  protected:
    Coord() = delete;
    explicit Coord(const size_t d): dimensions(d) {}

    template <class T> inline
    double _euclidean_distance(const std::vector<T>& v) const {
        double result = 0;
        for (const auto x: v) {result += x * x;}
        return std::sqrt(result);
    }

    std::vector<std::vector<int>> directions_;
    size_t max_neighbors_;
};

class Neumann final: public Coord {
  public:
    Neumann() = delete;
    explicit Neumann(const size_t d): Coord(d) {
        max_neighbors_ = 2 * d;
        directions_.reserve(max_neighbors_);
        std::vector<int> v(dimensions, 0);
        v.back() += 1;
        do {
            directions_.push_back(v);
        } while (std::next_permutation(v.begin(), v.end()));
        v.assign(dimensions, 0);
        v.front() -= 1;
        do {
            directions_.push_back(v);
        } while (std::next_permutation(v.begin(), v.end()));
    }
    ~Neumann() = default;
    //! Manhattan distance
    virtual size_t graph_distance(const std::vector<int>& v) const override {
        return std::accumulate(v.begin(), v.end(), 0, [](const int lhs, const int rhs) {
            return std::abs(lhs) + std::abs(rhs);
        });
    }
};

//! Neumann + diagonal cells
class Moore final: public Coord {
  public:
    Moore() = delete;
    explicit Moore(const size_t d): Coord(d) {
        max_neighbors_ = std::pow(3, dimensions) - 1;
        directions_.reserve(max_neighbors_);
        for (const int x: {-1, 0, 1}) {
            for (const int y: {-1, 0, 1}) {
                if (dimensions == 2) {
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
    }
    ~Moore() = default;
    //! Chebyshev/chessboard distance
    virtual size_t graph_distance(const std::vector<int>& v) const override {
        return std::abs(*std::max_element(v.begin(), v.end(), [](const int lhs, const int rhs) {
            return std::abs(lhs) < std::abs(rhs);
        }));
    }
};

class Hexagonal final: public Coord {
  public:
    Hexagonal() = delete;
    explicit Hexagonal(const size_t d): Coord(d) {
        max_neighbors_ = 6 * (d - 1);
        std::vector<int> v{-1, 0, 1};
        directions_.reserve(max_neighbors_);
        if (dimensions == 2) {
            do {
                directions_.push_back({v[0], v[1]});
            } while (std::next_permutation(v.begin(), v.end()));
        }
        else {
            do {
                directions_.push_back({v[0], v[1], 0});
            } while (std::next_permutation(v.begin(), v.end()));
            directions_.push_back({0, 0, -1});
            directions_.push_back({1, 0, -1});
            directions_.push_back({1, -1, -1});
            directions_.push_back({0, 0, 1});
            directions_.push_back({-1, 0, 1});
            directions_.push_back({-1, 1, 1});
        }
    }
    ~Hexagonal() = default;
    virtual size_t graph_distance(const std::vector<int>& v) const override {
        std::vector<size_t> absv;
        absv.reserve(v.size() * 2);
        for (auto x: v) {
            absv.push_back(std::abs(x));
        }
        absv.push_back(std::abs(v[0] + v[1]));
        if (v.size() == 3) {
            absv.push_back(std::abs(v[0] + v[2]));
        }
        return *std::max_element(absv.begin(), absv.end());
    }
    virtual double euclidean_distance(const std::vector<int>& v) const override {
        std::vector<double> true_pos(v.begin(), v.end());
        true_pos[1] += true_pos[0] * 0.5;
        true_pos[0] *= std::sqrt(3.0 / 4.0);
        if (v.size() > 2) {
            true_pos[0] += true_pos[2] / sqrt(3.0);
            true_pos[2] *= std::sqrt(2.0 / 3.0);
        }
        return _euclidean_distance(true_pos);
    }
    virtual double radius(size_t volume) const override {
        if (dimensions == 2) {
            return Coord::radius(volume) * std::sqrt(3.0 / 4.0);
        } else {
            return Coord::radius(volume) * std::sqrt(0.5);
        }
    }
    virtual std::vector<std::vector<int>> core() const override {
        std::vector<std::vector<int>> output = Neumann(dimensions).core();
        if (dimensions == 3) {
            output.resize(3);
            output.push_back({1, 0, -1});
        }
        return output;
    }
};



#endif /* COORD_HPP_ */
