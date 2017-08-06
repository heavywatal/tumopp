// -*- mode: c++; coding: utf-8 -*-
/*! @file tissue.hpp
    @brief Interface of Tissue class
*/
#pragma once
#ifndef TISSUE_HPP_
#define TISSUE_HPP_

#include "coord.hpp"
#include "cell.hpp"

#include <boost/functional/hash.hpp>

#include <cstdint>
#include <iostream>
#include <vector>
#include <valarray>
#include <unordered_set>
#include <memory>
#include <string>

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

namespace boost {
    namespace program_options {
        class options_description;
    }
}

namespace std {
  template <> struct hash<std::valarray<int>> {
    size_t operator() (const std::valarray<int>& v) const {
        return boost::hash_range(std::begin(v), std::end(v));
    }
  };
  template <> struct hash<std::shared_ptr<tumopp::Cell>> {
    size_t operator() (const std::shared_ptr<tumopp::Cell>& x) const {
        return std::hash<std::valarray<int>>()(x->coord());
    }
  };
}

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

namespace tumopp {

class equal_shptr_cell {
  public:
    bool operator() (const std::shared_ptr<tumopp::Cell>& lhs,
                     const std::shared_ptr<tumopp::Cell>& rhs) const {
        return (lhs->coord() == rhs->coord()).min();
    }
};

class Tissue {
  public:
    static unsigned int DIMENSIONS() {return DIMENSIONS_;}

    //! Constructor
    Tissue() = default;

    //! main function
    bool grow(const size_t max_size, const double plateau=0.0);

    std::ostream& write_segsites(std::ostream&, const std::vector<std::shared_ptr<Cell>>&) const;
    std::vector<std::shared_ptr<Cell>> sample_random(const size_t) const;
    std::vector<std::shared_ptr<Cell>> sample_section(const size_t) const;

    std::string specimens() const {return specimens_.str();}
    std::string snapshots() const {return snapshots_.str();}
    std::string drivers() const {return drivers_.str();}
    std::string pairwise_distance(const size_t npair) const;

    void write(std::ostream&) const;

    friend std::ostream& operator<< (std::ostream&, const Tissue&);

    size_t size() const {return tumor_.size();};
    double radius() const {return coord_func_->radius(tumor_.size());}
    double relative_pos(std::shared_ptr<Cell> x) const {
        return coord_func_->euclidean_distance(x->coord()) / radius();
    }

    //! Set coordinate function object
    template <class FuncObj>
    void set_coord() {coord_func_ = std::make_unique<FuncObj>(DIMENSIONS_);}
    const std::unique_ptr<Coord>& coord_func() const {return coord_func_;}

    //! Unit test
    static void unit_test();
    static boost::program_options::options_description opt_description();

  private:
    //! Dimensions: {1, 2, 3}
    static unsigned int DIMENSIONS_;

    //! Coordinate system
    static std::string COORDINATE_;

    //! {const, step, linear}
    static std::string LOCAL_DENSITY_EFFECT_;

    //! {random, mindrag, minstraight}
    static std::string DISPLACEMENT_PATH_;

    //! 0: flat, +: peripheral growth
    static double SIGMA_E_;

    //! initial population size
    static size_t INITIAL_SIZE_;

    //! a flag
    static size_t RECORDING_EARLY_GROWTH_;

    //! Time (tumor_.size()) to introduce a driver mutation
    static size_t MUTATION_TIMING_;

    void init();
    void init_coord();
    void init_insert_function();
    std::function<bool(const std::shared_ptr<Cell>&)> insert;

    //! Swap with a random neighbor
    void migrate(const std::shared_ptr<Cell>&);

    //! Emplace daughter cell and push other cells to the direction
    void push(std::shared_ptr<Cell> moving, const std::valarray<int>& direction);
    //! Push through the minimum drag path
    void push_minimum_drag(std::shared_ptr<Cell> moving);
    //! Try insert_adjacent() on every step in push()
    void stroll(std::shared_ptr<Cell> moving, const std::valarray<int>& direction);
    //! Insert x if any adjacent node is empty
    bool insert_adjacent(const std::shared_ptr<Cell>& x);
    //! Put new cell and return existing.
    bool swap_existing(std::shared_ptr<Cell>* x);
    //! Count steps to the nearest empty
    size_t steps_to_empty(std::valarray<int> current, const std::valarray<int>& direction) const;
    //! Direction to the nearest empty
    std::valarray<int> to_nearest_empty(const std::valarray<int>& current, const unsigned int search_max=26) const;
    //! Direction is selected with a probability proportional with 1/l
    std::valarray<int> roulette_direction(const std::valarray<int>& current) const;

    uint_fast8_t num_empty_neighbors(const std::valarray<int>&) const;
    double proportion_empty_neighbors(const std::valarray<int>& coord) const {
        double x = static_cast<double>(num_empty_neighbors(coord));
        return x /= static_cast<double>(coord_func_->directions().size());
    }

    void queue_push(const std::shared_ptr<Cell>&);

    double positional_value(const std::valarray<int>&) const;

    std::vector<size_t> generate_neutral_mutations() const;

    std::string header() const;
    void write(std::ostream&, const Cell&) const;

    /////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
    // Data member

    //! cells
    std::unordered_set<
        std::shared_ptr<Cell>,
        std::hash<std::shared_ptr<Cell>>,
        equal_shptr_cell> tumor_;

    //! event queue
    std::multimap<double, std::shared_ptr<Cell>> queue_;

    //! continuous time
    double time_ = 0.0;

    //! incremented when a new cell is born
    size_t id_tail_ = 0;

    //! Call set_coord() once
    std::unique_ptr<Coord> coord_func_;

    //! record all cells existed
    std::ostringstream specimens_;
    //! record state when size() < RECORDING_EARLY_GROWTH_
    std::ostringstream snapshots_;
    //! record driver mutations
    std::ostringstream drivers_;
};

} // namespace tumopp

#endif /* TISSUE_HPP_ */
