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
  //! hash function for cell coord
  template <> struct hash<std::valarray<int>> {
    //! hash function for cell coord
    size_t operator() (const std::valarray<int>& v) const {
        return boost::hash_range(std::begin(v), std::end(v));
    }
  };
  //! hash function for cell coord
  template <> struct hash<std::shared_ptr<tumopp::Cell>> {
    //! hash function for cell coord
    size_t operator() (const std::shared_ptr<tumopp::Cell>& x) const {
        return std::hash<std::valarray<int>>()(x->coord());
    }
  };
}

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

namespace tumopp {

//! Compare cell coord
class equal_shptr_cell {
  public:
    //! Compare cell coord
    bool operator() (const std::shared_ptr<tumopp::Cell>& lhs,
                     const std::shared_ptr<tumopp::Cell>& rhs) const {
        return (lhs->coord() == rhs->coord()).min();
    }
};

/*! a mass of cancer cells
*/
class Tissue {
  public:
    //! getter of DIMENSIONS_
    static unsigned int DIMENSIONS() {return DIMENSIONS_;}

    //! Constructor
    Tissue() = default;

    //! main function
    bool grow(const size_t max_size, const double plateau=0.0);

    //! Write ms-like binary sequence
    std::ostream& write_segsites(std::ostream&, const std::vector<std::shared_ptr<Cell>>&) const;
    //! sample cells
    std::vector<std::shared_ptr<Cell>> sample_random(const size_t) const;
    //! sample cells in a cross section
    std::vector<std::shared_ptr<Cell>> sample_section(const size_t) const;

    //! Stringify #specimens_
    std::string specimens() const {return specimens_.str();}
    //! Stringify #snapshots_
    std::string snapshots() const {return snapshots_.str();}
    //! Stringify #drivers_
    std::string drivers() const {return drivers_.str();}
    //! Make TSV of pairwise distance
    std::string pairwise_distance(const size_t npair) const;

    //! Write all cells
    void write(std::ostream&) const;

    friend std::ostream& operator<< (std::ostream&, const Tissue&);

    //! Shortcut of tumor_.size()
    size_t size() const {return tumor_.size();};
    //! Shortcut of coord.radius(tumor_.size())
    double radius() const {return coord_func_->radius(tumor_.size());}
    //! Relative distance from the origin
    double relative_pos(std::shared_ptr<Cell> x) const {
        return coord_func_->euclidean_distance(x->coord()) / radius();
    }

    //! Set #coord_func_ for testing
    template <class FuncObj>
    void init_coord_test() {coord_func_ = std::make_unique<FuncObj>(DIMENSIONS_);}
    //! getter of #coord_func_
    const std::unique_ptr<Coord>& coord_func() const {return coord_func_;}

    //! Unit test
    static void test();

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

    //! Initializer is separated from constructor for restarting
    void init();
    //! Set #coord_func_
    void init_coord();
    //! Check LOCAL_DENSITY_EFFECT_ and DISPLACEMENT_PATH_ and set function
    void init_insert_function();
    //! initialized in init_insert_function()
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

    //! Count adjacent empty sites
    uint_fast8_t num_empty_neighbors(const std::valarray<int>&) const;
    //! num_empty_neighbors() / directions().size()
    double proportion_empty_neighbors(const std::valarray<int>& coord) const {
        double x = static_cast<double>(num_empty_neighbors(coord));
        return x /= static_cast<double>(coord_func_->directions().size());
    }

    //! Push a cell to event queue
    void queue_push(const std::shared_ptr<Cell>&);

    //! Calculate positional value
    double positional_value(const std::valarray<int>&) const;

    //! Put random mutations on genealogy
    std::vector<size_t> generate_neutral_mutations() const;

    //! TSV header
    std::string header() const;
    //! TSV
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
    //! record state when size() < #RECORDING_EARLY_GROWTH_
    std::ostringstream snapshots_;
    //! record driver mutations
    std::ostringstream drivers_;
};

} // namespace tumopp

#endif /* TISSUE_HPP_ */
