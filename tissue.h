// -*- mode: c++; coding: utf-8 -*-
/*! @file tissue.h
    @brief Interface of Tissue class
*/
#pragma once
#ifndef TISSUE_H_
#define TISSUE_H_

#include <iostream>
#include <vector>
#include <unordered_set>
#include <memory>
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
  template <> struct hash<std::shared_ptr<Gland>> {
    size_t operator() (const std::shared_ptr<Gland>& x) const {
        return std::hash<std::vector<int>>()(x->coord());
    }
  };
}

class equal_shptr_gland {
  public:
    bool operator() (const std::shared_ptr<Gland>& lhs,
                     const std::shared_ptr<Gland>& rhs) {
        return lhs->coord() == rhs->coord();
    }
};

namespace boost {
    namespace program_options {
        class options_description;
    }
}


class Tissue {
  public:
    static size_t DIMENSIONS() {return DIMENSIONS_;}

    //! Constructor
    Tissue(): coord_func_(new Neumann(DIMENSIONS_)) {};

    //! Set coordinate function object
    template <class FuncObj>
    void set_coord() {coord_func_.reset(new FuncObj(DIMENSIONS_));}

    //! Mark first cells with mutation
    void stain();

    //!
    void grow(const size_t max_size);

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

    //! replication schedule: random, poisson, even
    static std::string SCHEDULE_;

    //! packing method: push, push_fill, walk_fill
    static std::string PACKING_;

    //! Emplace daughter gland and push other glands outward
    void push(const std::shared_ptr<Gland>& moving, std::vector<int> direction={});
    //! Push if not exists an empty neighbor
    void push_fill(const std::shared_ptr<Gland>& daughter, const std::vector<int>& direction={});
    //! Change direction before push_fill
    void walk_fill(const std::shared_ptr<Gland>& daughter);
    //! Fill inner layer first and move outward
    void push_layer(const std::shared_ptr<Gland>& daughter);
    //! Do not push and fill empty space
    bool fill_empty(const std::shared_ptr<Gland>& daughter);

    std::unique_ptr<Coord> coord_func_;

    std::unordered_set<std::shared_ptr<Gland>,
        std::hash<std::shared_ptr<Gland>>,
        equal_shptr_gland> tumor_;
    //! The coordinates of the past mutations
    std::vector<std::vector<int>> mutation_coords_;
    //! Timing of mutations (tumor size as proxy)
    std::vector<size_t> mutation_stages_;
    std::vector<std::string> evolution_history_;

    const std::string sep_ = "\t";
};

#endif /* TISSUE_H_ */
