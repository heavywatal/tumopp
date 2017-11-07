/*! @file coord.hpp
    @brief Coordinate system
*/
#pragma once
#ifndef COORD_HPP_
#define COORD_HPP_

#include <stdexcept>
#include <valarray>
#include <vector>
#include <random>

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
namespace tumopp {

//! constexpr pow function for integers
constexpr size_t ipow(size_t base, size_t exponent) noexcept {
  return exponent <= 0u ? 1u
       : exponent == 1u ? base
       : base * ipow(base, --exponent);
}

/*! @brief Base class of coordinate system
*/
class Coord {
  public:
    //! getter of #directions_
    const std::vector<std::valarray<int>>& directions() const {return directions_;}

    //! [0, 0, ...]
    std::valarray<int> origin() const {
        return std::valarray<int>(dimensions_);
    }
    //! Convert #directions_ to absolute coordinates
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
        return directions_[dist_direction_(rng)];
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
    double cross_section(size_t nodes) const;

    // virtual methods
    //! Graph distance
    virtual int graph_distance(const std::valarray<int>& v) const = 0;
    //! Euclidean distance
    virtual double euclidean_distance(const std::valarray<int>& v) const {
        return _euclidean_distance(v);
    }
    //! Estimate radius from volume
    virtual double radius(const size_t nodes) const;
    //! square or cube
    virtual std::vector<std::valarray<int>> core() const;
    //! sphere coordinates with inside-out direction
    std::vector<std::valarray<int>> sphere(const size_t n) const;
    //! Destructor
    virtual ~Coord() = default;

  protected:
    //! Default constructor is deleted
    Coord() = delete;
    //! Constructor: initialize and check #dimensions_
    explicit Coord(const unsigned int d): dimensions_(d) {
        if (d < 2U || 3U < d) {
            throw std::runtime_error("Invalid value for dimensions");
        }
    }

    //! Euclidean distance
    template <class T> inline
    double _euclidean_distance(const std::valarray<T>& v) const {
        return std::sqrt((v * v).sum());
    }

    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
    // Data member

    //! {2, 3}
    const unsigned int dimensions_;
    //! initialized in derived class constructor
    std::vector<std::valarray<int>> directions_;
    //! initialized in derived class constructor
    std::uniform_int_distribution<unsigned int> dist_direction_;
};

/*! @brief Derived class of Coord
*/
class Neumann final: public Coord {
  public:
    Neumann() = delete;
    //! Constructor
    explicit Neumann(const unsigned int d): Coord(d) {
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
    ~Neumann() = default;
    //! Manhattan distance
    int graph_distance(const std::valarray<int>& v) const override;
};

/*! @brief Derived class of Coord

    Neumann + diagonal cells
*/
class Moore final: public Coord {
  public:
    Moore() = delete;
    //! Constructor
    explicit Moore(const unsigned int d): Coord(d) {
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
    ~Moore() = default;
    //! Chebyshev/chessboard distance
    int graph_distance(const std::valarray<int>& v) const override;
};

/*! @brief Derived class of Coord
*/
class Hexagonal final: public Coord {
  public:
    Hexagonal() = delete;
    //! Constructor
    explicit Hexagonal(const unsigned int d): Coord(d) {
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
    ~Hexagonal() = default;
    int graph_distance(const std::valarray<int>& v) const override;
    double euclidean_distance(const std::valarray<int>& v) const override;
    double radius(const size_t nodes) const override;
    std::vector<std::valarray<int>> core() const override;
};

} // namespace tumopp

#endif /* COORD_HPP_ */
