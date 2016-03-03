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

#include "cell.h"
#include "coord.hpp"

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

namespace std {
  template <> struct hash<std::vector<int>> {
    size_t operator() (const std::vector<int>& v) const {
        return boost::hash_range(v.begin(), v.end());
    }
  };
  template <> struct hash<std::shared_ptr<Cell>> {
    size_t operator() (const std::shared_ptr<Cell>& x) const {
        return std::hash<std::vector<int>>()(x->coord());
    }
  };
}

class equal_shptr_cell {
  public:
    bool operator() (const std::shared_ptr<Cell>& lhs,
                     const std::shared_ptr<Cell>& rhs) const {
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
    Tissue() {
        if (COORDINATE_ == "neumann") {coord_func_ = std::make_unique<Neumann>(DIMENSIONS_);}
        else if (COORDINATE_ == "moore") {coord_func_ = std::make_unique<Moore>(DIMENSIONS_);}
        else if (COORDINATE_ == "hex") {coord_func_ = std::make_unique<Hexagonal>(DIMENSIONS_);}
    }

    //! Set coordinate function object
    template <class FuncObj>
    void set_coord() {coord_func_ = std::make_unique<FuncObj>(DIMENSIONS_);}
    const std::unique_ptr<Coord>& coord_func() const {return coord_func_;}

    //! main function
    void grow(const size_t max_size);

    std::ostream& write_segsites(std::ostream&, const std::vector<std::shared_ptr<Cell>>&) const;
    std::vector<std::shared_ptr<Cell>> sample_random(const size_t) const;
    std::vector<std::shared_ptr<Cell>> sample_if(std::function<bool(const std::vector<int>&)>) const;

    //! Return tumor state as TSV string
    std::string snapshot() const;
    //! Return mutation_coords_ as TSV string
    std::string mutation_history() const;

    std::string snapshot_header() const;
    std::string evolution_history() const;

    friend std::ostream& operator<< (std::ostream&, const Tissue&);

    size_t size() const {return tumor_.size();};
    double radius() const {return coord_func_->radius(tumor_.size());}
    double relative_pos(std::shared_ptr<Cell> x) const {
        return coord_func_->euclidean_distance(x->coord()) / radius();
    }

    //! Unit test
    static void unit_test();
    static boost::program_options::options_description& opt_description();

  private:
    //! Dimensions: {1, 2, 3}
    static size_t DIMENSIONS_;

    //! Coordinate system
    static std::string COORDINATE_;

    //! packing method: push, push_fill, walk_fill
    static std::string PACKING_;

    //! 0: flat, +: peripheral growth
    static double GLOBAL_ENV_COEF_;

    bool insert(const std::shared_ptr<Cell>&);

    //! Swap with a random neighbor
    void migrate(const std::shared_ptr<Cell>&);

    //! Emplace daughter cell and push other cells outward
    void push(std::shared_ptr<Cell> moving, const std::vector<int>& direction);
    void pushn_everytime(std::shared_ptr<Cell> moving);
    //! Fill empty neighbor or push to the direction
    void fill_push(std::shared_ptr<Cell> moving, const std::vector<int>& direction);
    //! Insert x if it has an empty neighbor
    bool fill_empty(const std::shared_ptr<Cell>& x);
    //! Try to insert x to a random neighbor
    bool insert_neighbor(const std::shared_ptr<Cell>& daughter);

    //! Put new cell and return existing.
    bool swap_existing(std::shared_ptr<Cell>* x);
    std::vector<std::vector<int>> empty_neighbors(const std::vector<int>&) const;
    size_t steps_to_empty(std::vector<int> current, const std::vector<int>& direction) const;
    std::vector<int> to_nearest_empty(const std::vector<int>& current, size_t search_max=26) const;

    void queue_push(double t, const std::shared_ptr<Cell>&);

    double positional_value(const std::vector<int>&) const;

    std::unique_ptr<Coord> coord_func_;

    std::unordered_set<std::shared_ptr<Cell>,
        std::hash<std::shared_ptr<Cell>>,
        equal_shptr_cell> tumor_;

    //! all existed cells
    std::vector<std::shared_ptr<Cell>> stock_;

    //! event queue
    std::multimap<double, std::shared_ptr<Cell>> queue_;

    //! The coordinates of the past mutations
    std::vector<std::vector<int>> mutation_coords_;

    //! Timing of mutations (tumor size as proxy)
    std::vector<size_t> mutation_stages_;
    std::vector<std::string> evolution_history_;

    const std::string sep_ = "\t";
};

#endif /* TISSUE_H_ */
