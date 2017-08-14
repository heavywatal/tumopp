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

namespace boost {
    namespace program_options {
        class options_description;
    }
}

namespace tumopp {

/*! @brief Population of Cell
*/
class Tissue {
  public:
    //! Constructor
    Tissue() = default;

    //! main function
    bool grow(const size_t max_size, const double plateau=0.0);

    //! getter of #DIMENSIONS_
    static unsigned int DIMENSIONS() {return DIMENSIONS_;}
    //! Shortcut of tumor_.size()
    size_t size() const {return tumor_.size();};

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

    //! Unit test
    static void test();
    //! Shortcut of coord.radius(tumor_.size())
    double radius() const {return coord_func_->radius(tumor_.size());}
    //! getter of #coord_func_
    const std::unique_ptr<Coord>& coord_func() const {return coord_func_;}
    //! Set #coord_func_ for testing
    template <class FuncObj>
    void init_coord_test() {coord_func_ = std::make_unique<FuncObj>(DIMENSIONS_);}

    static boost::program_options::options_description opt_description();

  private:
    //! Initializer is separated from constructor for restarting
    void init();
    //! Set #coord_func_
    void init_coord();
    //! Set #insert function according to #LOCAL_DENSITY_EFFECT_ and #DISPLACEMENT_PATH_
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

    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
    // Function object for tumor_

    //! Hashing function object for cell coord
    struct hash_valarray_int {
        //! hash function
        size_t operator() (const std::valarray<int>& v) const {
            return boost::hash_range(std::begin(v), std::end(v));
        }
    };

    //! Hashing function object for shptr<Cell>
    struct hash_shptr_cell {
        //! hash function
        size_t operator() (const std::shared_ptr<tumopp::Cell>& x) const {
            return hash_valarray_int()(x->coord());
        }
    };

    //! Equal function object for shptr<Cell>
    struct equal_shptr_cell {
        //! Compare cell coord
        bool operator() (const std::shared_ptr<tumopp::Cell>& lhs,
                         const std::shared_ptr<tumopp::Cell>& rhs) const {
            return (lhs->coord() == rhs->coord()).min();
        }
    };

    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
    // Data member

    //! Dimensions: {2, 3}
    static unsigned int DIMENSIONS_;
    //! Coordinate/neighborhood system {neumann, moore, hex}
    static std::string COORDINATE_;
    //! E2 {const, step, linear}
    static std::string LOCAL_DENSITY_EFFECT_;
    //! Push method {1: random, 2: roulette, 3: mindrag, 4: minstraight, 5: stroll}
    static std::string DISPLACEMENT_PATH_;
    //! 0: flat, +: peripheral growth
    static double SIGMA_E_;
    //! initial population size
    static size_t INITIAL_SIZE_;
    //! a flag
    static size_t RECORDING_EARLY_GROWTH_;
    //! Time (tumor_.size()) to introduce a driver mutation
    static size_t MUTATION_TIMING_;

    //! cells
    std::unordered_set<
        std::shared_ptr<Cell>,
        hash_shptr_cell,
        equal_shptr_cell> tumor_;
    //! incremented when a new cell is born
    size_t id_tail_ = 0;

    //! event queue
    std::multimap<double, std::shared_ptr<Cell>> queue_;
    //! continuous time
    double time_ = 0.0;

    //! initialized in init_coord() or init_coord_test()
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
