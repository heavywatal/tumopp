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
      const EventRates& init_event_rates=EventRates{},
      const CellParams& cell_params=CellParams{});

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

    //! Put random mutations on genealogy
    std::vector<unsigned> generate_neutral_mutations(double mu, bool has_at_least_1_mutation_per_division) const;
    //! Write ms-like binary sequence
    std::ostream& write_segsites(std::ostream&, const std::vector<std::shared_ptr<Cell>>&, const std::vector<unsigned>&) const;
    //! sample cells around the given cell
    std::vector<std::shared_ptr<Cell>> sample_bulk(const std::shared_ptr<Cell>&, size_t) const;
    //! sample cells uniformly
    std::vector<std::shared_ptr<Cell>> sample_medoids(size_t) const;
    //! sample cells randomly
    std::vector<std::shared_ptr<Cell>> sample_random(size_t) const;
    //! sample cells in a cross section
    std::vector<std::shared_ptr<Cell>> sample_section(size_t) const;

    //! Stringify #specimens_
    std::string specimens();
    //! Stringify #snapshots_
    std::string snapshots() const {
        return "time\t" + header() + snapshots_.str();
    }
    //! Stringify #drivers_
    std::string drivers() const {
        return "id\ttype\tcoef\n" + drivers_.str();
    }
    //! Make TSV of pairwise distance
    std::string pairwise_distance(size_t npair) const;

    //! TSV header
    std::string header() const;
    //! TSV
    void write(std::ostream&, const Cell&) const;
    //! Write all cells to #snapshots_ with #time_
    void write_snapshot();
    friend std::ostream& operator<< (std::ostream&, const Tissue&);

    //! @cond
    bool snapshots_empty() const {return snapshots_.str().empty();};
    bool drivers_empty() const {return drivers_.str().empty();}
    //! @endcond

    //! Shortcut of tumor_.size()
    size_t size() const noexcept {return tumor_.size();}
    //! Shortcut of tumor_.begin()
    auto begin() noexcept {return tumor_.begin();}
    //! Shortcut of tumor_.end()
    auto end() noexcept {return tumor_.end();}
    //! getter of #time_
    double time() const noexcept {return time_;}
    //! Shortcut of coord_func_.dimensions()
    unsigned dimensions() const noexcept {return coord_func_->dimensions();}

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

    //! TODO: Calculate positional value
    double positional_value(const std::valarray<int>&) const {return 1.0;}

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
    //! counter of periodic snapshots
    unsigned i_snapshot_ = 1u;
    //! initialized in init_coord() or init_coord_test()
    std::unique_ptr<Coord> coord_func_;

    //! record all cells existed
    std::vector<std::shared_ptr<Cell>> specimens_;
    //! record snapshots
    std::ostringstream snapshots_;
    //! record driver mutations
    std::ostringstream drivers_;
};

} // namespace tumopp

#endif /* TISSUE_HPP_ */
