// -*- mode: c++; coding: utf-8 -*-
/*! @file tissue.h
    @brief Interface of Tissue class
*/
#pragma once
#ifndef TISSUE_H_
#define TISSUE_H_
#include <cmath>
#include <cassert>
#include <iostream>
#include <vector>
#include <map>
#include <unordered_map>
#include <string>

#include <boost/functional/hash.hpp>

#include "gland.h"

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


namespace std {
  template <> struct hash<std::vector<int>> {
    size_t operator() (const std::vector<int>& v) const {
        return boost::hash_range(v.begin(), v.end());
    }
  };
}

namespace boost {
    namespace program_options {
        class options_description;
    }
}

class Orthogonal;
class Lattice;
class Hex;

class Tissue {
    typedef Hex coord_t;
  public:
    static size_t DIMENSIONS() {return DIMENSIONS_;}

    //! Constructor
    Tissue(const std::vector<int>& origin=std::vector<int>(DIMENSIONS_)) {
        emplace(origin, Gland());
        init_regularly();
    }

    //! Initiate tumor with neighboring glands regularly
    void init_regularly();
    //! Initiate tumor with neighboring glands randomly
    void init_randomly();
    //! Mark first cells with mutation
    void stain();

    //! Big bang model
    void grow_random(const size_t max_size);
    //! Another extreme
    void grow_even(const size_t max_size);
    //! in-between
    void grow_poisson(const size_t max_size);

    //! Return tumor state as TSV string
    std::string snapshot() const;
    //! Return mutation_coords_ as TSV string
    std::string mutation_history() const;

    std::string snapshot_header() const;
    std::string evolution_history() const;

    friend std::ostream& operator<< (std::ostream&, const Tissue&);

    //! Unit test
    static void unit_test();
    static boost::program_options::options_description& opt_description();

  private:
    //! Dimensions: {1, 2, 3}
    static size_t DIMENSIONS_;

    //! Emplace daughter gland at the specified coord
    void emplace(const std::vector<int>& current_coord, Gland&& daughter);
    //! Emplace daughter gland and push other glands outward
    void push(Gland&& daughter, std::vector<int>* current_coord,
              const std::vector<int>& direction);
    //! Push if not exists an empty neighbor
    void push_fill(Gland&& daughter, const std::vector<int>& current_coord,
                   const std::vector<int>& direction={});
    //! Change direction before push_fill
    void walk_fill(Gland&& daughter, const std::vector<int>& current_coord);

    //! key: coords, value: gland
    std::unordered_map<std::vector<int>, Gland> tumor_;
    //! The coordinates of the existing glands
    std::vector<std::vector<int>> coords_;
    //! The coordinates of the past mutations
    std::vector<std::vector<int>> mutation_coords_;
    //! Timing of mutations (tumor size as proxy)
    std::vector<size_t> mutation_stages_;
    std::vector<std::string> evolution_history_;

    const std::string sep_ = "\t";
};

#endif /* TISSUE_H_ */
