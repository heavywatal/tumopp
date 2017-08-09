// -*- mode: c++; coding: utf-8 -*-
/*! @file coord.hpp
    @brief Coordinate system
*/
#pragma once
#ifndef COORD_HPP_
#define COORD_HPP_

#include <boost/math/constants/constants.hpp>

#include <cmath>
#include <cstdlib> // std::abs for int
#include <stdexcept>
#include <valarray>
#include <vector>
#include <bitset>
#include <numeric>
#include <random>

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
namespace tumopp {

//! Mathematical constant
constexpr double PI = boost::math::constants::pi<double>();

/*! Base class of coordinate system
*/
class Coord {
  public:
    //! getter of #directions_
    const std::vector<std::valarray<int>>& directions() const {return directions_;}

    //! [0, 0, ...]
    std::valarray<int> origin() const {
        return std::valarray<int>(dimensions);
    }
    //! List neighboring sites
    std::vector<std::valarray<int>> neighbors(const std::valarray<int>& v) const {
        std::vector<std::valarray<int>> output = directions_;
        for (auto& d: output) {
            d += v;
        }
        return output;
    }
    //! Choose a random neighbor
    template <class RNG> inline
    std::valarray<int> random_direction(RNG& rng) {
        return directions_[uniform_int_dist_(rng)];
    }
    //! Choose a random neighbor of the specified site
    template <class RNG> inline
    std::valarray<int> random_neighbor(const std::valarray<int>& v, RNG& rng) {
        return v + random_direction(rng);
    }
    //! Direction that maximize the distance from the origin
    std::valarray<int> outward(const std::valarray<int>& v) const {
        const auto candidates = neighbors(v);
        return *std::max_element(candidates.begin(), candidates.end(),
                                 [this](const std::valarray<int>& lhs, const std::valarray<int>& rhs) {
            return euclidean_distance(lhs) < euclidean_distance(rhs);
        });
    }
    //! Area of cross section
    double cross_section(size_t vol) const {
        return std::pow((9.0 / 16.0) * PI * (vol *= vol), 1.0 / 3.0);
    }

    // virtual methods
    //! Graph distance
    virtual int graph_distance(const std::valarray<int>& v) const = 0;
    //! Euclidean distance
    virtual double euclidean_distance(const std::valarray<int>& v) const {
        return _euclidean_distance(v);
    }
    //! Estimate radius from volume
    virtual double radius(const size_t points) const {
        double x = points;
        x /= PI;
        if (dimensions == 2U) {
            // S = pi r^2
            return std::sqrt(x);
        } else {
            // V = 4/3 pi r^3
            return std::pow(x *= (3.0 / 4.0), 1.0 / 3.0);
        }
    }
    //! square or cube
    virtual std::vector<std::valarray<int>> core() const {
        const size_t n = std::pow(2, dimensions);
        std::vector<std::valarray<int>> output;
        output.reserve(n);
        for (size_t i=0; i<n; ++i) {
            std::bitset<3> bs(i);
            std::valarray<int> v(dimensions);
            for (size_t j=0; j<dimensions; ++j) {
                v[j] = static_cast<int>(bs[j]);
            }
            output.push_back(v);
        }
        return output;
    }
    //! sphere coordinates with inside-out direction
    std::vector<std::valarray<int>> sphere(const size_t n) const {
        std::vector<std::valarray<int>> output;
        if (dimensions == 2U) {
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
    //! Destructor
    virtual ~Coord() = default;

    //! two or three
    const unsigned int dimensions;

  protected:
    //! Default constructor is deleted
    Coord() = delete;
    //! Initialize with the number of dimensions
    explicit Coord(const unsigned int d): dimensions(d) {
        if (d < 2U || 3U < d) {
            throw std::runtime_error("Invalid value for dimensions");
        }
    }

    //! Euclidean distance
    template <class T> inline
    double _euclidean_distance(const std::valarray<T>& v) const {
        return std::sqrt((v * v).sum());
    }

    //! initialized in derived class constructor
    std::vector<std::valarray<int>> directions_;
    //! initialized in derived class constructor
    std::uniform_int_distribution<ptrdiff_t> uniform_int_dist_;
};

/*! Derived class of coordinate system
*/
class Neumann final: public Coord {
  public:
    Neumann() = delete;
    //! Constructor
    explicit Neumann(const size_t d): Coord(d) {
        directions_.reserve(2U * dimensions);
        std::valarray<int> v(dimensions);
        v[v.size() - 1] += 1;
        do {
            directions_.push_back(v);
        } while (std::next_permutation(std::begin(v), std::end(v)));
        v = 0;
        v[0] -= 1;
        do {
            directions_.push_back(v);
        } while (std::next_permutation(std::begin(v), std::end(v)));
        uniform_int_dist_ = std::uniform_int_distribution<ptrdiff_t>(0, directions_.size() - 1);
    }
    ~Neumann() = default;
    //! Manhattan distance
    virtual int graph_distance(const std::valarray<int>& v) const override {
        const std::valarray<int> absv = std::abs(v);
        return absv.sum();
    }
};

/*! Derived class of coordinate system

    Neumann + diagonal cells
*/
class Moore final: public Coord {
  public:
    Moore() = delete;
    //! Constructor
    explicit Moore(const size_t d): Coord(d) {
        directions_.reserve(std::pow(3, dimensions) - 1);
        for (const int x: {-1, 0, 1}) {
            for (const int y: {-1, 0, 1}) {
                if (dimensions == 2U) {
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
        uniform_int_dist_ = std::uniform_int_distribution<ptrdiff_t>(0, directions_.size() - 1);
    }
    ~Moore() = default;
    //! Chebyshev/chessboard distance
    virtual int graph_distance(const std::valarray<int>& v) const override {
        return std::abs(v).max();
    }
};

/*! Derived class of coordinate system
*/
class Hexagonal final: public Coord {
  public:
    Hexagonal() = delete;
    //! Constructor
    explicit Hexagonal(const size_t d): Coord(d) {
        std::valarray<int> v{-1, 0, 1};
        directions_.reserve(6 * (dimensions - 1));
        if (dimensions == 2U) {
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
        uniform_int_dist_ = std::uniform_int_distribution<ptrdiff_t>(0, directions_.size() - 1);
    }
    ~Hexagonal() = default;
    virtual int graph_distance(const std::valarray<int>& v) const override {
        int d = std::max(std::abs(v).max(), std::abs(v[0] + v[1]));
        if (dimensions > 2U) {
            return std::max(d, std::abs(v[0] + v[2]));
        }
        return d;
    }
    virtual double euclidean_distance(const std::valarray<int>& v) const override {
        std::valarray<double> true_pos(dimensions);
        true_pos[0] += static_cast<double>(v[0]);
        true_pos[1] += static_cast<double>(v[1]);
        true_pos[1] += true_pos[0] * 0.5;
        true_pos[0] *= std::sqrt(3.0 / 4.0);
        if (dimensions > 2U) {
            true_pos[2] += static_cast<double>(v[2]);
            true_pos[0] += true_pos[2] / sqrt(3.0);
            true_pos[2] *= std::sqrt(2.0 / 3.0);
        }
        return _euclidean_distance(true_pos);
    }
    virtual double radius(const size_t volume) const override {
        if (dimensions == 2U) {
            return Coord::radius(volume * std::sqrt(3.0 / 4.0));
        } else {
            return Coord::radius(volume * std::sqrt(0.5));
        }
    }
    virtual std::vector<std::valarray<int>> core() const override {
        std::vector<std::valarray<int>> output = Neumann(dimensions).core();
        if (dimensions == 3U) {
            output.resize(3);
            output.push_back({1, 0, -1});
        }
        return output;
    }
};

} // namespace tumopp

#endif /* COORD_HPP_ */
