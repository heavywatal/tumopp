/*! @file tissue.hpp
    @brief Interface of Tissue class
*/
#pragma once
#ifndef TUMOPP_TISSUE_HPP_
#define TUMOPP_TISSUE_HPP_

#include "coord.hpp"
#include "cell.hpp"
#include "random.hpp"

#include <cstdint>
#include <string>
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
      ptrdiff_t initial_size=1,
      unsigned dimensions=3u,
      const std::string& coordinate="moore",
      const std::string& local_density_effect="const",
      const std::string& displacement_path="random",
      const EventRates& init_event_rates=EventRates{},
      uint32_t seed=std::random_device{}(),
      bool verbose=false,
      bool enable_benchmark=false);
    ~Tissue();

    //! main function
    bool grow(
      ptrdiff_t max_size,
      double max_time=100.0,
      double snapshot_interval=0.0,
      ptrdiff_t recording_early_growth=0,
      ptrdiff_t mutation_timing=0);

    //! Simulate turnover with the increased death_rate
    void plateau(double time);
    //! Simulate medical treatment with the increased death_prob
    void treatment(double death_prob, ptrdiff_t num_resistant_cells = 3);

    //! Append current state to #benchmark_
    void benchmark_append(const ptrdiff_t size);
    //! Write #extant_cells_ and their ancestors
    std::string str_history() const;
    //! Write #snapshots_
    std::string str_snapshots() const;
    //! Write #drivers_
    std::string str_drivers() const;
    //! View #benchmark_
    std::string_view str_benchmark() const { return benchmark_; }
    //! Write all cells to buffer
    std::string format() const;

    //! @cond
    bool has_snapshots() const {return !snapshots_.empty();};
    bool has_drivers() const {return !drivers_.empty();}
    //! @endcond

    //! @name Getter functions
    //@{
    //! Get the number of extant cells
    size_t size() const noexcept {return extant_cells_.size();}
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
    ptrdiff_t steps_to_empty(const coord_t& current, const coord_t& direction) const;
    //! Direction to the nearest empty
    const coord_t& to_nearest_empty(const coord_t& current) const;
    //! Direction is selected with a probability proportional with 1/l
    coord_t roulette_direction(const coord_t& current) const;

    //! Count adjacent empty sites
    int_fast8_t num_empty_neighbors(const coord_t&) const;
    //! TODO: Calculate positional value
    double positional_value(const coord_t&) const {return 1.0;}

    //! Push a cell to event #queue_
    void queue_push(const std::shared_ptr<Cell>&, bool surrounded=false);
    //! Put a cell to #cemetery_
    void entomb(const std::shared_ptr<Cell>&);
    //! Write all cells to #snapshots_ with #time_
    void snapshots_append();

    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
    // Function object for extant_cells_

    //! Hashing function object for shared_pointer<Cell>
    struct hash_ptr_cell {
        //! hash function
        size_t operator() (const std::shared_ptr<tumopp::Cell>& x) const noexcept {
            return hash(x->coord());
        }
    };

    //! Equal function object for shared_pointer<Cell>
    struct equal_ptr_cell {
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
        hash_ptr_cell,
        equal_ptr_cell> extant_cells_{};
    //! incremented when a new cell is born
    int id_tail_{0};

    //! event queue
    std::multimap<double, std::shared_ptr<Cell>> queue_{};
    //! continuous time
    double time_{0.0};
    //! initialized in init_coord() or init_coord_test()
    std::unique_ptr<Coord> coord_func_{nullptr};

    //! record dead cells
    std::string cemetery_{};
    //! record snapshots
    std::string snapshots_{};
    //! record driver mutations
    std::string drivers_{};
    //! id of recorded cells
    mutable std::unordered_set<int> recorded_{};
    //! record resource usage
    std::string benchmark_{};
    //! random number generator
    std::unique_ptr<urbg_t> engine_{nullptr};
    //! print debug info
    bool verbose_{false};
};

inline std::string format_as(const Tissue& x) { return x.format(); }

} // namespace tumopp

#endif // TUMOPP_TISSUE_HPP_
