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
#include <array>
#include <unordered_set>
#include <map>
#include <memory>
#include <functional>

namespace tumopp {

/*! @brief Population of Cell
*/
class Tissue {
  public:
    //! Constructor
    Tissue(
      size_t initial_size=1u,
      unsigned dimensions=3u,
      const std::string& coordinate="moore",
      const std::string& local_density_effect="const",
      const std::string& displacement_path="random",
      const EventRates& init_event_rates=EventRates{});

    //! main function
    bool grow(
      size_t max_size,
      double max_time=100.0,
      double snapshot_interval=std::numeric_limits<double>::infinity(),
      size_t recording_early_growth=0u,
      size_t mutation_timing=std::numeric_limits<size_t>::max());

    //! Simulate turnover with the increased death_rate
    void plateau(double time);
    //! Simulate medical treatment with the increased death_prob
    void treatment(double death_prob, size_t num_resistant_cells = 3u);

    //! Stringify #extant_cells_ and their ancestors
    std::stringstream history() const;
    //! Stringify #snapshots_
    std::stringstream snapshots() const;
    //! Stringify #drivers_
    std::stringstream drivers() const;
    friend std::ostream& operator<< (std::ostream&, const Tissue&);

    //! @cond
    bool has_snapshots() const {return snapshots_.rdbuf()->in_avail();};
    bool has_drivers() const {return drivers_.rdbuf()->in_avail();}
    //! @endcond

    //! @name Getter functions
    size_t size() const noexcept {return extant_cells_.size();}
    auto begin() noexcept {return extant_cells_.begin();}
    auto end() noexcept {return extant_cells_.end();}
    unsigned dimensions() const noexcept {return coord_func_->dimensions();}
    //@}

  private:
    //! Set #coord_func_
    void init_coord(unsigned dimensions, const std::string& coordinate);
    //! Set #insert function
    void init_insert_function(const std::string& local_density_effect, const std::string& displacement_path);
    //! initialized in init_insert_function()
    std::function<bool(const std::shared_ptr<Cell>&)> insert;

    //! Swap with a random neighbor
    void migrate(const std::shared_ptr<Cell>&);
    //! Emplace daughter cell and push other cells to the direction
    void push(std::shared_ptr<Cell> moving, const coord_t& direction);
    //! Push through the minimum drag path
    void push_minimum_drag(std::shared_ptr<Cell> moving);
    //! Try insert_adjacent() on every step in push()
    void stroll(std::shared_ptr<Cell> moving, const coord_t& direction);
    //! Insert x if any adjacent node is empty
    bool insert_adjacent(const std::shared_ptr<Cell>& x);
    //! Put new cell and return existing.
    bool swap_existing(std::shared_ptr<Cell>* x);
    //! Count steps to the nearest empty
    size_t steps_to_empty(coord_t current, const coord_t& direction) const;
    //! Direction to the nearest empty
    coord_t to_nearest_empty(const coord_t& current, const unsigned int search_max=26) const;
    //! Direction is selected with a probability proportional with 1/l
    coord_t roulette_direction(const coord_t& current) const;

    //! Count adjacent empty sites
    uint_fast8_t num_empty_neighbors(const coord_t&) const;
    //! num_empty_neighbors() / directions().size()
    double proportion_empty_neighbors(const coord_t& coord) const {
        double x = static_cast<double>(num_empty_neighbors(coord));
        return x /= static_cast<double>(coord_func_->directions().size());
    }
    //! TODO: Calculate positional value
    double positional_value(const coord_t&) const {return 1.0;}

    //! Push a cell to event #queue_
    void queue_push(const std::shared_ptr<Cell>&, bool surrounded=false);
    //! Write all cells to #snapshots_ with #time_
    void snapshots_append();

    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
    // Function object for extant_cells_

    //! Hashing function object for shptr<Cell>
    struct hash_shptr_cell {
        //! hash function
        size_t operator() (const std::shared_ptr<tumopp::Cell>& x) const noexcept {
            return hash(x->coord());
        }
    };

    //! Equal function object for shptr<Cell>
    struct equal_shptr_cell {
        //! Compare cell coord
        bool operator() (const std::shared_ptr<tumopp::Cell>& lhs,
                         const std::shared_ptr<tumopp::Cell>& rhs) const noexcept {
            return lhs->coord() == rhs->coord();
        }
    };

    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
    // Data member

    //! cells
    std::unordered_set<
        std::shared_ptr<Cell>,
        hash_shptr_cell,
        equal_shptr_cell> extant_cells_;
    //! incremented when a new cell is born
    unsigned id_tail_ = 0;

    //! event queue
    std::multimap<double, std::shared_ptr<Cell>> queue_;
    //! continuous time
    double time_ = 0.0;
    //! counter of periodic snapshots
    unsigned i_snapshot_ = 1u;
    //! initialized in init_coord() or init_coord_test()
    std::unique_ptr<Coord> coord_func_;

    //! record snapshots
    std::stringstream snapshots_;
    //! record driver mutations
    std::stringstream drivers_;
};

} // namespace tumopp

#endif /* TISSUE_HPP_ */
