/*! @file tissue.cpp
    @brief Implementation of Tissue class
*/
#include "tissue.hpp"
#include "random.hpp"

#include <wtl/random.hpp>
#include <wtl/iostr.hpp>
#include <wtl/debug.hpp>
#include <wtl/math.hpp>
#include <wtl/numeric.hpp>
#include <wtl/algorithm.hpp>
#include <wtl/genetic.hpp>

namespace tumopp {

Tissue::Tissue(
  const size_t initial_size,
  const unsigned dimensions,
  const std::string& coordinate,
  const std::string& local_density_effect,
  const std::string& displacement_path,
  const EventRates& init_event_rates) {HERE;
    snapshots_.precision(std::cout.precision());
    drivers_.precision(std::cout.precision());
    init_coord(dimensions, coordinate);
    init_insert_function(local_density_effect, displacement_path);
    const auto initial_coords = coord_func_->sphere(initial_size);
    const auto origin = std::make_shared<Cell>(
      initial_coords[0], ++id_tail_,
      std::make_shared<EventRates>(init_event_rates)
    );
    extant_cells_.insert(origin);
    while (extant_cells_.size() < initial_size) {
        for (const auto& mother: extant_cells_) {
            const auto daughter = std::make_shared<Cell>(*mother);
            const auto ancestor = std::make_shared<Cell>(*mother);
            ancestor->set_time_of_death(0.0);
            mother->set_time_of_birth(0.0, ++id_tail_, ancestor);
            daughter->set_time_of_birth(0.0, ++id_tail_, ancestor);
            daughter->set_coord(initial_coords[extant_cells_.size()]);
            extant_cells_.insert(daughter);
            if (extant_cells_.size() >= initial_size) break;
        }
    }
    for (const auto& cell: extant_cells_) queue_push(cell);
}

void Tissue::init_coord(const unsigned dimensions, const std::string& coordinate) {HERE;
    std::unordered_map<std::string, std::unique_ptr<Coord>> swtch;
    swtch["neumann"] = std::make_unique<Neumann>(dimensions);
    swtch["moore"] = std::make_unique<Moore>(dimensions);
    swtch["hex"] = std::make_unique<Hexagonal>(dimensions);
    try {
        coord_func_ = std::move(swtch.at(coordinate));
    } catch (std::exception& e) {
        std::ostringstream oss;
        oss << std::endl << FILE_LINE_PRETTY
            << "\nInvalid value for -C (" << coordinate << "); choose from "
            << wtl::keys(swtch);
        throw std::runtime_error(oss.str());
    }
}

bool Tissue::grow(const size_t max_size, const double max_time,
                  const double snapshot_interval,
                  size_t recording_early_growth,
                  size_t mutation_timing) {HERE;
    if (recording_early_growth > 0u) {snapshots_append();}
    bool success = false;
    size_t i = 0;
    double time_snapshot = i_snapshot_ * snapshot_interval;
    while (true) {
        if ((++i % 1000U) == 0U) {DCERR("\r" << extant_cells_.size());}
        auto it = queue_.begin();
        time_ = it->first;
        if (time_ > max_time || extant_cells_.size() >= max_size) {
            success = true; // maybe not; but want to exit with record
            break;
        }
        if (time_snapshot > 0.0 && time_ > time_snapshot) {
            snapshots_append();
            time_snapshot = ++i_snapshot_ * snapshot_interval;
        }
        const auto mother = std::move(it->second);
        queue_.erase(it);
        if (mother->next_event() == Event::birth) {
            const auto daughter = std::make_shared<Cell>(*mother);
            if (insert(daughter)) {
                const auto ancestor = std::make_shared<Cell>(*mother);
                ancestor->set_time_of_death(time_);
                mother->set_time_of_birth(time_, ++id_tail_, ancestor);
                daughter->differentiate();
                daughter->set_time_of_birth(time_, ++id_tail_, ancestor);
                drivers_ << mother->mutate();
                drivers_ << daughter->mutate();
                if (extant_cells_.size() == mutation_timing) {
                    mutation_timing = 0u; // once
                    drivers_ << daughter->force_mutate();
                }
                queue_push(mother);
                queue_push(daughter);
            } else {
                queue_push(mother, true);
                continue;  // skip write()
            }
        } else if (mother->next_event() == Event::death) {
            extant_cells_.erase(mother);
            if (extant_cells_.empty()) break;
        } else {
            migrate(mother);
            queue_push(mother);
        }
        if (extant_cells_.size() < recording_early_growth) {
            snapshots_append();
        } else {
            recording_early_growth = 0u;  // prevent restart by cell death
        }
    }
    DCERR("\r" << extant_cells_.size() << std::endl);
    return success;
}

void Tissue::plateau(const double time) {HERE;
    queue_.clear();
    for (auto& p: extant_cells_) {
        p->increase_death_rate();
        queue_push(p);
    }
    grow(std::numeric_limits<size_t>::max(), time_ + time);
}

void Tissue::treatment(const double death_prob, const size_t num_resistant_cells) {HERE;
    const size_t original_size = extant_cells_.size();
    std::vector<std::shared_ptr<Cell>> cells;
    cells.reserve(original_size);
    for (const auto& p: queue_) { // for reproducibility
        cells.emplace_back(p.second);
    }
    std::shuffle(cells.begin(), cells.end(), engine64());
    for (size_t i=0; i<original_size; ++i) {
        const auto& p = cells[i];
        if (i >= num_resistant_cells) {
            p->set_cycle_dependent_death(death_prob);
        }
    }
    const size_t margin = 10u * num_resistant_cells + 10u;
    grow(original_size + margin, std::numeric_limits<double>::max());
}

void Tissue::queue_push(const std::shared_ptr<Cell>& x, const bool surrounded) {
    double dt = x->delta_time(time_, positional_value(x->coord()), surrounded);
    queue_.emplace_hint(queue_.end(), dt += time_, x);
}

void Tissue::init_insert_function(const std::string& local_density_effect, const std::string& displacement_path) {
    using func_t = std::function<bool(const std::shared_ptr<Cell>&)>;
    using map_sf = std::unordered_map<std::string, func_t>;
    std::unordered_map<std::string, map_sf> swtch;

    swtch["const"].emplace("random", [this](const std::shared_ptr<Cell>& daughter) {
        push(daughter, coord_func_->random_direction(engine64()));
        return true;
    });
    swtch["const"].emplace("mindrag", [this](const std::shared_ptr<Cell>& daughter) {
        push_minimum_drag(daughter);
        return true;
    });
    swtch["const"].emplace("minstraight", [this](const std::shared_ptr<Cell>& daughter) {
        push(daughter, to_nearest_empty(daughter->coord()));
        return true;
    });
    swtch["const"].emplace("roulette", [this](const std::shared_ptr<Cell>& daughter) {
        push(daughter, roulette_direction(daughter->coord()));
        return true;
    });
    swtch["const"].emplace("stroll", [this](const std::shared_ptr<Cell>& daughter) {
        stroll(daughter, coord_func_->random_direction(engine64()));
        return true;
    });
    swtch["step"].emplace("random", [this](const std::shared_ptr<Cell>& daughter) {
        if (num_empty_neighbors(daughter->coord()) == 0U) {return false;}
        push(daughter, coord_func_->random_direction(engine64()));
        return true;
    });
    swtch["step"].emplace("mindrag", [this](const std::shared_ptr<Cell>& daughter) {
        return insert_adjacent(daughter);
    });
    swtch["linear"].emplace("random", [this](const std::shared_ptr<Cell>& daughter) {
        const double prob = proportion_empty_neighbors(daughter->coord());
        if (wtl::generate_canonical(engine64()) < prob) {
            push(daughter, coord_func_->random_direction(engine64()));
            return true;
        }
        return false;
    });
    swtch["linear"].emplace("mindrag", [this](const std::shared_ptr<Cell>& daughter) {
        daughter->set_coord(coord_func_->random_neighbor(daughter->coord(), engine64()));
        return extant_cells_.insert(daughter).second;
    });
    try {
        insert = swtch.at(local_density_effect).at(displacement_path);
    } catch (std::exception& e) {
        std::ostringstream oss;
        oss << std::endl << FILE_LINE_PRETTY
            << "\nInvalid value for -L/-P ("
            << local_density_effect << "/" << displacement_path
            << "); choose from";
        for (const auto& p: swtch) {
            oss << "\n -L" << p.first << " -P " << wtl::keys(p.second);
        }
        throw std::runtime_error(oss.str());
    }
}

void Tissue::push(std::shared_ptr<Cell> moving, const coord_t& direction) {
    do {
        moving->set_coord(moving->coord() + direction);
    } while (swap_existing(&moving));
}

void Tissue::push_minimum_drag(std::shared_ptr<Cell> moving) {
    do {
        moving->set_coord(moving->coord() + to_nearest_empty(moving->coord()));
    } while (swap_existing(&moving));
}

void Tissue::stroll(std::shared_ptr<Cell> moving, const coord_t& direction) {
    while (!insert_adjacent(moving)) {
        moving->set_coord(moving->coord() + direction);
        swap_existing(&moving);
    }
}

bool Tissue::insert_adjacent(const std::shared_ptr<Cell>& moving) {
    const auto present_coord = moving->coord();
    auto neighbors = coord_func_->neighbors(present_coord);
    std::shuffle(neighbors.begin(), neighbors.end(), engine64());
    for (auto& x: neighbors) {
        moving->set_coord(x);
        if (extant_cells_.insert(moving).second) {
            return true;
        }
    }
    moving->set_coord(present_coord);
    return false;
}

bool Tissue::swap_existing(std::shared_ptr<Cell>* x) {
    // The cell must not be in extant_cells_.
    auto result = extant_cells_.insert(*x);
    if (result.second) {
        return false;
    } else {
        std::shared_ptr<Cell> existing = std::move(*result.first);
        extant_cells_.erase(result.first);
        extant_cells_.insert(std::move(*x));
        x->swap(existing);
        return true;
    }
}

void Tissue::migrate(const std::shared_ptr<Cell>& moving) {
    extant_cells_.erase(moving);
    const auto orig_pos = moving->coord();
    moving->set_coord(coord_func_->random_neighbor(moving->coord(), engine64()));
    auto result = extant_cells_.insert(moving);
    if (!result.second) {
        std::shared_ptr<Cell> existing = std::move(*result.first);
        extant_cells_.erase(result.first);
        extant_cells_.insert(std::move(moving));
        existing->set_coord(orig_pos);
        extant_cells_.insert(existing);
    }
}

size_t Tissue::steps_to_empty(coord_t current, const coord_t& direction) const {
    size_t steps = 0;
    const auto key = std::make_shared<Cell>();
    do {
        key->set_coord(current += direction);
        ++steps;
    } while (extant_cells_.find(key) != extant_cells_.end());
    return steps;
}

coord_t Tissue::to_nearest_empty(const coord_t& current, const unsigned search_max) const {
    size_t least_steps = std::numeric_limits<size_t>::max();
    coord_t best_direction{};
    auto directions = coord_func_->directions();
    std::shuffle(directions.begin(), directions.end(), engine64());
    if (search_max < directions.size()) directions.resize(search_max);
    for (const auto& d: directions) {
        auto n = steps_to_empty(current, d);
        if (n < least_steps) {
            least_steps = n;
            best_direction = d;
        }
    }
    return best_direction;
}

coord_t Tissue::roulette_direction(const coord_t& current) const {
    auto directions = coord_func_->directions();
    std::shuffle(directions.begin(), directions.end(), engine64());
    std::vector<double> roulette;
    for (const auto& d: directions) {
        const auto l = steps_to_empty(current, d);
        if (l == 0U) {return d;}
        roulette.push_back(1.0 / l);
    }
    return directions[wtl::roulette_select(roulette, engine64())];
}

uint_fast8_t Tissue::num_empty_neighbors(const coord_t& coord) const {
    uint_fast8_t cnt = 0;
    std::shared_ptr<Cell> nb = std::make_shared<Cell>();
    for (const auto& d: coord_func_->directions()) {
        nb->set_coord(coord + d);
        if (extant_cells_.find(nb) == extant_cells_.end()) {++cnt;}
    }
    return cnt;
}

std::stringstream Tissue::history() const {
    std::stringstream ss;
    ss.precision(std::cout.precision());
    ss << Cell::header() << "\n";
    std::unordered_set<unsigned> done;
    for (const auto& p: extant_cells_) {
        p->traceback(ss, &done);
    }
    return ss;
}

std::stringstream Tissue::snapshots() const {
    std::stringstream ss;
    ss << "time\t" << Cell::header() << "\n" << snapshots_.rdbuf();
    return ss;
}

std::stringstream Tissue::drivers() const {
    std::stringstream ss;
    ss << "id\ttype\tcoef\n" << drivers_.rdbuf();
    return ss;
}

void Tissue::snapshots_append() {
    for (const auto& p: extant_cells_) {
        snapshots_ << time_ << "\t" << *p << "\n";
    }
}

//! Stream operator for debug print
std::ostream& operator<< (std::ostream& ost, const Tissue& tissue) {
    for (const auto& p: tissue.extant_cells_) {
        ost << *p << "\n";
    }
    return ost;
}

} // namespace tumopp
