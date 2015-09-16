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
#include "coord.hpp"

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////


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


class Tissue {
  public:
    static size_t DIMENSIONS() {return DIMENSIONS_;}

    //! Constructor
    Tissue(const std::vector<int>& origin=std::vector<int>(DIMENSIONS_));

    //! Set coordinate function object
    template <class FuncObj>
    void set_coord() {coord_func_.reset(new FuncObj());}

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
              const std::vector<int>& direction={});
    //! Push if not exists an empty neighbor
    void push_fill(Gland&& daughter, const std::vector<int>& current_coord,
                   const std::vector<int>& direction={});
    //! Change direction before push_fill
    void walk_fill(Gland&& daughter, const std::vector<int>& current_coord);

    std::unique_ptr<Coord> coord_func_;

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
