/*! @file tissue.hpp
    @brief Interface of Tissue class
*/
#pragma once
#ifndef TISSUE_HPP_
#define TISSUE_HPP_

#include "coord.hpp"
#include "cell.hpp"

#include <cstdint>
#include <sstream>
#include <string>
#include <vector>
#include <valarray>
#include <unordered_set>
#include <map>
#include <memory>
#include <functional>

namespace boost {namespace program_options {class options_description;}}

namespace tumopp {

/*! @brief Population of Cell
*/
class Tissue {
  public:
    //! Constructor
    Tissue() = default;

    //! main function
    bool grow(size_t max_size, double plateau=0.0);

    //! getter of #DIMENSIONS_
    static unsigned int DIMENSIONS() noexcept {return DIMENSIONS_;}
    //! Shortcut of tumor_.size()
    size_t size() const noexcept {return tumor_.size();}
    //! Shortcut of tumor_.begin()
    auto begin() noexcept {return tumor_.begin();}
    //! Shortcut of tumor_.end()
    auto end() noexcept {return tumor_.end();}

    //! Write ms-like binary sequence
    std::ostream& write_segsites(std::ostream&, const std::vector<std::shared_ptr<Cell>>&) const;
    //! sample cells around the given cell
    std::vector<std::shared_ptr<Cell>> sample_bulk(const std::shared_ptr<Cell>&, size_t) const;
    //! sample cells uniformly
    std::vector<std::shared_ptr<Cell>> sample_medoids(size_t) const;
    //! sample cells randomly
    std::vector<std::shared_ptr<Cell>> sample_random(size_t) const;
    //! sample cells in a cross section
    std::vector<std::shared_ptr<Cell>> sample_section(size_t) const;

    //! Stringify #specimens_
    std::string specimens() const {return specimens_.str();}
    //! Stringify #snapshots_
    std::string snapshots() const {return snapshots_.str();}
    //! Stringify #drivers_
    std::string drivers() const {return drivers_.str();}
    //! Make TSV of pairwise distance
    std::string pairwise_distance(size_t npair) const;

    //! TSV header
    std::string header() const;
    //! TSV
    void write(std::ostream&, const Cell&) const;
    //! Write all cells
    void write_snapshot();
    friend std::ostream& operator<< (std::ostream&, const Tissue&);

    //! @cond
    bool has_snapshots() const;
    bool has_drivers() const;
    //! @endcond

    //! Shortcut of coord.radius(tumor_.size())
    double radius() const {return coord_func_->radius(tumor_.size());}
    //! getter of #coord_func_
    const std::unique_ptr<Coord>& coord_func() const noexcept {return coord_func_;}

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
    std::vector<unsigned> generate_neutral_mutations() const;

    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
    // Function object for tumor_

    //! Hashing function object for shptr<Cell>
    struct hash_shptr_cell {
        //! hash function
        size_t operator() (const std::shared_ptr<tumopp::Cell>& x) const noexcept {
            return Coord::hash(x->coord());
        }
    };

    //! Equal function object for shptr<Cell>
    struct equal_shptr_cell {
        //! Compare cell coord
        bool operator() (const std::shared_ptr<tumopp::Cell>& lhs,
                         const std::shared_ptr<tumopp::Cell>& rhs) const noexcept {
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
    unsigned id_tail_ = 0;

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
