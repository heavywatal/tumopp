/*! @file coord.hpp
    @brief Coordinate system
*/
#pragma once
#ifndef COORD_HPP_
#define COORD_HPP_

#include <array>
#include <vector>
#include <random>

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
namespace tumopp {

//! Maximum dimensions
constexpr unsigned MAX_DIM = 3u;

//! Alias of coordinate type
using coord_t = std::array<int, MAX_DIM>;

//! hash coordinates for std::unordered_* container
inline size_t hash(const coord_t& v) {
    // no collision as long as each of v is within +-1e6
    constexpr coord_t::value_type mask21 = 0x001fffff;
    size_t seed = v[0u] & mask21;
    seed |= static_cast<size_t>(v[1u] & mask21) << 21;
    seed |= static_cast<size_t>(v[2u] & mask21) << 42;
    return seed;
}


//! @name Arithmetic operators for std::array
//@{
template <class T> inline
std::array<T, MAX_DIM>& operator+=(std::array<T, MAX_DIM>& lhs, const std::array<T, MAX_DIM>& rhs) {
    for (unsigned i = 0u; i < MAX_DIM; ++i) {lhs[i] += rhs[i];}
    return lhs;
}

template <class T> inline
std::array<T, MAX_DIM>& operator-=(std::array<T, MAX_DIM>& lhs, const std::array<T, MAX_DIM>& rhs) {
    for (unsigned i = 0u; i < MAX_DIM; ++i) {lhs[i] -= rhs[i];}
    return lhs;
}

template <class T> inline
std::array<T, MAX_DIM>& operator*=(std::array<T, MAX_DIM>& lhs, const std::array<T, MAX_DIM>& rhs) {
    for (unsigned i = 0u; i < MAX_DIM; ++i) {lhs[i] *= rhs[i];}
    return lhs;
}

template <class T> inline
std::array<T, MAX_DIM> operator+(const std::array<T, MAX_DIM>& lhs, const std::array<T, MAX_DIM>& rhs) {
    std::array<T, MAX_DIM> v(lhs);
    return v += rhs;
}

template <class T> inline
std::array<T, MAX_DIM> operator-(const std::array<T, MAX_DIM>& lhs, const std::array<T, MAX_DIM>& rhs) {
    std::array<T, MAX_DIM> v(lhs);
    return v -= rhs;
}

template <class T> inline
std::array<T, MAX_DIM> operator*(const std::array<T, MAX_DIM>& lhs, const std::array<T, MAX_DIM>& rhs) {
    std::array<T, MAX_DIM> v(lhs);
    return v *= rhs;
}
//@}

/*! @brief Base class of coordinate system
*/
class Coord {
  public:
    //! Convert #directions_ to absolute coordinates
    std::vector<coord_t> neighbors(const coord_t& v) const {
        std::vector<coord_t> output = directions_;
        for (auto& d: output) {
            d += v;
        }
        return output;
    }
    //! Choose a random neighbor
    template <class URBG> inline
    coord_t random_direction(URBG& engine) {
        return directions_[dist_direction_(engine)];
    }
    //! Choose a random neighbor of the specified site
    template <class URBG> inline
    coord_t random_neighbor(const coord_t& v, URBG& engine) {
        return v + random_direction(engine);
    }
    //! Direction that maximize the distance from the origin
    coord_t outward(const coord_t& v) const;

    // virtual methods
    //! Convert coordinates into continuous scale
    virtual std::array<double, MAX_DIM> continuous(const coord_t& v) const;
    //! Graph distance
    virtual int graph_distance(const coord_t& v) const = 0;
    //! Euclidean distance
    virtual double euclidean_distance(const coord_t& v) const;
    //! square or cube
    virtual std::vector<coord_t> core() const;
    //! sphere coordinates with inside-out direction
    std::vector<coord_t> sphere(size_t n) const;
    //! Destructor
    virtual ~Coord() = default;

    //! getter of #directions_
    const std::vector<coord_t>& directions() const noexcept {return directions_;}
    //! getter of #dimensions_
    unsigned dimensions() const noexcept {return dimensions_;}

  protected:
    //! Default constructor is deleted
    Coord() = delete;
    //! Constructor: initialize and check #dimensions_
    explicit Coord(unsigned d);

    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
    // Data member

    //! {2, 3}
    const unsigned dimensions_;
    //! initialized in derived class constructor
    std::vector<coord_t> directions_;
    //! initialized in derived class constructor
    std::uniform_int_distribution<unsigned> dist_direction_;
};

/*! @brief Derived class of Coord
*/
class Neumann final: public Coord {
  public:
    Neumann() = delete;
    //! Constructor
    explicit Neumann(const unsigned d);
    ~Neumann() = default;
    //! Manhattan distance
    int graph_distance(const coord_t& v) const override;
};

/*! @brief Derived class of Coord

    Neumann + diagonal cells
*/
class Moore final: public Coord {
  public:
    Moore() = delete;
    //! Constructor
    explicit Moore(const unsigned d);
    ~Moore() = default;
    //! Chebyshev/chessboard distance
    int graph_distance(const coord_t& v) const override;
};

/*! @brief Derived class of Coord
*/
class Hexagonal final: public Coord {
  public:
    Hexagonal() = delete;
    //! Constructor
    explicit Hexagonal(const unsigned d);
    ~Hexagonal() = default;
    std::array<double, MAX_DIM> continuous(const coord_t& v) const override;
    int graph_distance(const coord_t& v) const override;
    double euclidean_distance(const coord_t& v) const override;
    std::vector<coord_t> core() const override;
};

} // namespace tumopp

#endif /* COORD_HPP_ */
