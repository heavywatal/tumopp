/*! @file coord.hpp
    @brief Coordinate system
*/
#pragma once
#ifndef TUMOPP_COORD_HPP_
#define TUMOPP_COORD_HPP_

#include <array>
#include <vector>
#include <random>

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
namespace tumopp {

//! Maximum dimensions
constexpr int MAX_DIM = 3;

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


//! @cond
//! @name Arithmetic operators for std::array
//@{
template <class T, class Func, std::size_t ...I> inline
std::array<T, sizeof...(I)>& array_arith_assign(
  std::array<T, sizeof...(I)>& lhs, const std::array<T, sizeof...(I)>& rhs,
  Func&& func, std::index_sequence<I...>) {
    (func(lhs[I], rhs[I]), ...);
    return lhs;
}

template <class T, class Func, std::size_t ...I> inline
std::array<T, sizeof...(I)>& array_val_arith_assign(
  std::array<T, sizeof...(I)>& lhs, const T rhs,
  Func&& func, std::index_sequence<I...>) {
    (func(lhs[I], rhs), ...);
    return lhs;
}

template <class T, std::size_t N> inline
std::array<T, N>& operator+=(std::array<T, N>& lhs, const std::array<T, N>& rhs) {
    return array_arith_assign(lhs, rhs, [](T& l, T r) { l += r; }, std::make_index_sequence<N>{});
}

template <class T, std::size_t N> inline
std::array<T, N>& operator-=(std::array<T, N>& lhs, const std::array<T, N>& rhs) {
    return array_arith_assign(lhs, rhs, [](T& l, T r) { l -= r; }, std::make_index_sequence<N>{});
}

template <class T, std::size_t N> inline
std::array<T, N>& operator*=(std::array<T, N>& lhs, const std::array<T, N>& rhs) {
    return array_arith_assign(lhs, rhs, [](T& l, T r) { l *= r; }, std::make_index_sequence<N>{});
}

template <class T, std::size_t N> inline
std::array<T, N>& operator*=(std::array<T, N>& lhs, T rhs) {
    return array_val_arith_assign(lhs, rhs, [](T& l, T r) { l *= r; }, std::make_index_sequence<N>{});
}

template <class T, std::size_t N> inline
std::array<T, N> operator+(const std::array<T, N>& lhs, const std::array<T, N>& rhs) {
    std::array<T, N> v(lhs);
    return v += rhs;
}

template <class T, std::size_t N> inline
std::array<T, N> operator-(const std::array<T, N>& lhs, const std::array<T, N>& rhs) {
    std::array<T, N> v(lhs);
    return v -= rhs;
}

template <class T, std::size_t N> inline
std::array<T, N> operator*(const std::array<T, N>& lhs, const std::array<T, N>& rhs) {
    std::array<T, N> v(lhs);
    return v *= rhs;
}

template <class T, std::size_t N> inline
std::array<T, N> operator*(const std::array<T, N>& lhs, T rhs) {
    std::array<T, N> v(lhs);
    return v *= rhs;
}

//@}
//! @endcond

/*! @brief Base class of coordinate system
*/
class Coord {
  public:
    //! Choose a random neighbor
    template <class URBG> inline
    coord_t random_direction(URBG& engine) {
        return directions_[dist_direction_(engine)];
    }

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
    std::vector<coord_t> sphere(ptrdiff_t n) const;
    //! Destructor
    virtual ~Coord() = default;

    //! getter of #directions_
    const std::vector<coord_t>& directions() const noexcept {return directions_;}
    //! getter of #dimensions_
    int dimensions() const noexcept {return dimensions_;}

  protected:
    //! Default constructor is deleted
    Coord() = delete;
    //! Constructor: initialize and check #dimensions_
    explicit Coord(int d);

    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
    // Data member

    //! {1, 2, 3}
    const int dimensions_{};
    //! initialized in derived class constructor
    std::vector<coord_t> directions_{};
    //! initialized in derived class constructor
    std::uniform_int_distribution<unsigned> dist_direction_{};
};

/*! @brief Derived class of Coord
*/
class Neumann final: public Coord {
  public:
    Neumann() = delete;
    //! Constructor
    explicit Neumann(const int d);
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
    explicit Moore(const int d);
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
    explicit Hexagonal(const int d);
    ~Hexagonal() = default;
    std::array<double, MAX_DIM> continuous(const coord_t& v) const override;
    int graph_distance(const coord_t& v) const override;
    double euclidean_distance(const coord_t& v) const override;
    std::vector<coord_t> core() const override;
};

} // namespace tumopp

#endif // TUMOPP_COORD_HPP_
