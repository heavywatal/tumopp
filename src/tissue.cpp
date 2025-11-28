/*! @file tissue.cpp
    @brief Implementation of Tissue class
*/
#include "tissue.hpp"

#include <fmt/format.h>
#include <fmt/ranges.h>
#include <wtl/random.hpp>
#include <wtl/numeric.hpp>
#include <wtl/algorithm.hpp>
#include <wtl/resource.hpp>
#include <wtl/signed.hpp>

#include <cstdio>

namespace tumopp {

Tissue::Tissue(
  const ptrdiff_t initial_size,
  const int dimensions,
  const std::string& coordinate,
  const std::string& local_density_effect,
  const std::string& displacement_path,
  const EventRates& init_event_rates,
  const uint32_t seed,
  const bool verbose,
  const bool enable_benchmark):
  engine_(std::make_unique<urbg_t>(seed)),
  verbose_(verbose) {
    if (enable_benchmark) {
        benchmark_append(0);
    }
    init_coord(dimensions, coordinate);
    init_insert_function(local_density_effect, displacement_path);
    const auto initial_coords = coord_func_->sphere(initial_size);
    const auto origin = std::make_shared<Cell>(
      initial_coords[0], ++id_tail_,
      std::make_shared<EventRates>(init_event_rates)
    );
    extant_cells_.insert(origin);
    while (wtl::ssize(extant_cells_) < initial_size) {
        for (const auto& mother: extant_cells_) {
            const auto daughter = std::make_shared<Cell>(*mother);
            const auto ancestor = std::make_shared<Cell>(*mother);
            ancestor->set_time_of_death(0.0);
            mother->set_time_of_birth(0.0, ++id_tail_, ancestor);
            daughter->set_time_of_birth(0.0, ++id_tail_, ancestor);
            daughter->set_coord(initial_coords[extant_cells_.size()]);
            extant_cells_.insert(daughter);
            if (wtl::ssize(extant_cells_) >= initial_size) break;
        }
    }
    for (const auto& cell: extant_cells_) queue_push(cell);
}

Tissue::~Tissue() = default;

void Tissue::init_coord(const int dimensions, const std::string& coordinate) {
    std::unordered_map<std::string, std::unique_ptr<Coord>> swtch;
    swtch["neumann"] = std::make_unique<Neumann>(dimensions);
    swtch["moore"] = std::make_unique<Moore>(dimensions);
    swtch["hex"] = std::make_unique<Hexagonal>(dimensions);
    try {
        coord_func_ = std::move(swtch.at(coordinate));
    } catch (std::exception& e) {
        const auto msg = fmt::format(
          "\n{}:{}:{}\nInvalid value for -C ({}); choose from {}\n",
          __FILE__, __LINE__, __PRETTY_FUNCTION__,
          coordinate, wtl::keys(swtch)
        );
        throw std::runtime_error(msg);
    }
}

bool Tissue::grow(const ptrdiff_t max_size, const double max_time,
                  const double snapshot_interval,
                  ptrdiff_t recording_early_growth,
                  ptrdiff_t mutation_timing) {
    if (recording_early_growth > 0) {snapshots_append();}
    bool success = false;
    double time_snapshot = snapshot_interval;
    constexpr ptrdiff_t progress_interval{1 << 12};
    while (true) {
        auto it = queue_.begin();
        time_ = it->first;
        if (time_ > max_time || wtl::ssize(extant_cells_) >= max_size) {
            success = true; // maybe not; but want to exit with record
            break;
        }
        if (time_snapshot > 0.0 && time_ > time_snapshot) {
            snapshots_append();
            time_snapshot = time_ + snapshot_interval;
        }
        const auto mother = std::move(it->second);
        queue_.erase(it);
        if (mother->next_event() == Event::birth) {
            const auto daughter = std::make_shared<Cell>(*mother);
            if (insert(daughter)) {
                const auto ancestor = std::make_shared<Cell>(*mother);
                ancestor->set_time_of_death(time_);
                mother->set_time_of_birth(time_, ++id_tail_, ancestor);
                daughter->differentiate(*engine_);
                daughter->set_time_of_birth(time_, ++id_tail_, ancestor);
                mother->mutate(drivers_, *engine_);
                daughter->mutate(drivers_, *engine_);
                if (wtl::ssize(extant_cells_) == mutation_timing) {
                    mutation_timing = 0; // once
                    daughter->force_mutate(drivers_, *engine_);
                }
                queue_push(mother);
                queue_push(daughter);
                const auto size = wtl::ssize(extant_cells_);
                if ((size % progress_interval) == 0) {
                    if (verbose_) fmt::print(stderr, "\r{}", size);
                    if (!benchmark_.empty()) benchmark_append(size);
                }
            } else {
                queue_push(mother, true);
                continue;  // skip write()
            }
        } else if (mother->next_event() == Event::death) {
            entomb(mother);
            if (extant_cells_.empty()) break;
        } else {
            migrate(mother);
            queue_push(mother);
        }
        if (wtl::ssize(extant_cells_) < recording_early_growth) {
            snapshots_append();
        } else {
            recording_early_growth = 0;  // prevent restart by cell death
        }
    }
    if (verbose_) fmt::println(stderr, "\r{}", extant_cells_.size());
    return success;
}

void Tissue::plateau(const double time) {
    queue_.clear();
    for (auto& p: extant_cells_) {
        p->increase_death_rate();
        queue_push(p);
    }
    grow(std::numeric_limits<ptrdiff_t>::max(), time_ + time);
}

void Tissue::treatment(const double death_prob, const ptrdiff_t num_resistant_cells) {
    const ptrdiff_t original_size = wtl::ssize(extant_cells_);
    std::vector<std::shared_ptr<Cell>> cells;
    wtl::reserve(cells, original_size);
    for (const auto& [_time, sh_ptr]: queue_) { // for reproducibility
        cells.emplace_back(sh_ptr);
    }
    std::shuffle(cells.begin(), cells.end(), *engine_);
    for (ptrdiff_t i=0; i<original_size; ++i) {
        const auto& p = cells[wtl::cast_u(i)];
        if (i >= num_resistant_cells) {
            p->set_cycle_dependent_death(*engine_, death_prob);
        }
    }
}

void Tissue::queue_push(const std::shared_ptr<Cell>& x, const bool surrounded) {
    double dt = x->delta_time(*engine_, time_, positional_value(x->coord()), surrounded);
    queue_.emplace_hint(queue_.end(), dt += time_, x);
}

void Tissue::init_insert_function(const std::string& local_density_effect, const std::string& displacement_path) {
    using func_t = std::function<bool(const std::shared_ptr<Cell>&)>;
    using map_sf = std::unordered_map<std::string, func_t>;
    std::unordered_map<std::string, map_sf> swtch;

    swtch["const"].emplace("random", [this](const std::shared_ptr<Cell>& daughter) {
        push(daughter, coord_func_->random_direction(*engine_));
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
        stroll(daughter, coord_func_->random_direction(*engine_));
        return true;
    });
    swtch["step"].emplace("random", [this](const std::shared_ptr<Cell>& daughter) {
        if (num_empty_neighbors(daughter->coord()) == 0) {return false;}
        push(daughter, coord_func_->random_direction(*engine_));
        return true;
    });
    swtch["step"].emplace("mindrag", [this](const std::shared_ptr<Cell>& daughter) {
        return insert_adjacent(daughter);
    });
    swtch["linear"].emplace("random", [this](const std::shared_ptr<Cell>& daughter) {
        const auto x = num_empty_neighbors(daughter->coord());
        if (x > 0) {
            double prob = x;
            prob /= static_cast<double>(coord_func_->directions().size());
            if (wtl::generate_canonical(*engine_) < prob) {
                push(daughter, coord_func_->random_direction(*engine_));
                return true;
            }
        }
        return false;
    });
    swtch["linear"].emplace("mindrag", [this](const std::shared_ptr<Cell>& daughter) {
        daughter->add_coord(coord_func_->random_direction(*engine_));
        return extant_cells_.insert(daughter).second;
    });
    try {
        insert = swtch.at(local_density_effect).at(displacement_path);
    } catch (std::exception& e) {
        std::string msg;
        fmt::format_to(std::back_inserter(msg),
          "\n{}:{}:{}\nInvalid value for -L/-P ({}/{}); choose from",
          __FILE__, __LINE__, __PRETTY_FUNCTION__,
          local_density_effect, displacement_path
        );
        for (const auto& [key, map]: swtch) {
            fmt::format_to(std::back_inserter(msg), "\n -L{} -P {}", key, wtl::keys(map));
        }
        throw std::runtime_error(msg);
    }
}

void Tissue::push(std::shared_ptr<Cell> moving, const coord_t& direction) {
    do {
        moving->add_coord(direction);
    } while (swap_existing(&moving));
}

void Tissue::push_minimum_drag(std::shared_ptr<Cell> moving) {
    do {
        moving->add_coord(to_nearest_empty(moving->coord()));
    } while (swap_existing(&moving));
}

void Tissue::stroll(std::shared_ptr<Cell> moving, const coord_t& direction) {
    while (!insert_adjacent(moving)) {
        moving->add_coord(direction);
        swap_existing(&moving);
    }
}

bool Tissue::insert_adjacent(const std::shared_ptr<Cell>& moving) {
    const auto present_coord = moving->coord();
    const auto& directions = coord_func_->directions();
    thread_local auto indices = wtl::seq_len<unsigned>(directions.size());
    std::shuffle(indices.begin(), indices.end(), *engine_);
    for (const auto i: indices) {
        moving->add_coord(directions[i]);
        if (extant_cells_.insert(moving).second) {
            return true;
        }
        moving->set_coord(present_coord);
    }
    return false;
}

bool Tissue::swap_existing(std::shared_ptr<Cell>* x) {
    auto [it, success] = extant_cells_.insert(*x);
    if (success) {
        return false;
    } else {
        std::shared_ptr<Cell> existing(*it);
        extant_cells_.insert(extant_cells_.erase(it), std::move(*x));
        *x = std::move(existing);
        return true;
    }
}

void Tissue::migrate(const std::shared_ptr<Cell>& migrant) {
    extant_cells_.erase(migrant);
    auto orig_pos = migrant->coord();
    migrant->add_coord(coord_func_->random_direction(*engine_));
    auto [it, success] = extant_cells_.insert(migrant);
    if (!success) {
        std::shared_ptr<Cell> existing(*it);
        extant_cells_.insert(extant_cells_.erase(it), migrant);
        existing->set_coord(std::move(orig_pos));
        extant_cells_.insert(std::move(existing));
    }
}

ptrdiff_t Tissue::steps_to_empty(const coord_t& current, const coord_t& direction) const {
    thread_local const auto key = std::make_shared<Cell>();
    const auto& end = extant_cells_.end();
    key->set_coord(current);
    ptrdiff_t steps = 0;
    do {
        key->add_coord(direction);
        ++steps;
    } while (extant_cells_.find(key) != end);
    return steps;
}

const coord_t& Tissue::to_nearest_empty(const coord_t& current) const {
    thread_local const auto key = std::make_shared<Cell>();
    const auto& end = extant_cells_.end();
    const auto& directions = coord_func_->directions();
    thread_local auto indices = wtl::seq_len<unsigned>(directions.size());
    std::shuffle(indices.begin(), indices.end(), *engine_);
    for (int radius = 1; true; ++radius) {
        for (const auto i: indices) {
            key->set_coord(current + directions[i] * radius);
            if (extant_cells_.find(key) == end) {
                return directions[i];
            }
        }
    }
}

coord_t Tissue::roulette_direction(const coord_t& current) const {
    auto directions = coord_func_->directions();
    std::shuffle(directions.begin(), directions.end(), *engine_);
    std::vector<double> roulette;
    for (const auto& d: directions) {
        const auto l = steps_to_empty(current, d);
        if (l == 0) {return d;}
        roulette.push_back(1.0 / static_cast<double>(l));
    }
    std::discrete_distribution<unsigned> discrete(roulette.begin(), roulette.end());
    return directions[discrete(*engine_)];
}

int_fast8_t Tissue::num_empty_neighbors(const coord_t& coord) const {
    thread_local auto nb = std::make_shared<Cell>();
    const auto& end = extant_cells_.end();
    int_fast8_t cnt = 0;
    for (const auto& d: coord_func_->directions()) {
        nb->set_coord(coord + d);
        if (extant_cells_.find(nb) == end) {++cnt;}
    }
    return cnt;
}

void Tissue::entomb(const std::shared_ptr<Cell>& dead) {
    dead->set_time_of_death(time_);
    dead->traceback(cemetery_, recorded_);
    extant_cells_.erase(dead);
}

std::string Tissue::str_history() const {
    std::string buffer = fmt::format("{}\n{}", Cell::header(), cemetery_);
    for (const auto& p: extant_cells_) {
        p->traceback(buffer, recorded_);
    }
    return buffer;
}

std::string Tissue::str_snapshots() const {
    return fmt::format("time\t{}\n{}", Cell::header(), snapshots_);
}

std::string Tissue::str_drivers() const {
    return fmt::format("id\ttype\tcoef\n{}", drivers_);
}

void Tissue::benchmark_append(const ptrdiff_t size) {
    const auto ru = wtl::get_rusage<std::milli, std::kilo>();
    fmt::format_to(std::back_inserter(benchmark_), "{}\t{}\n", size, fmt::join(ru, "\t"));
}

void Tissue::snapshots_append() {
    for (const auto& p: extant_cells_) {
        fmt::format_to(std::back_inserter(snapshots_), "{}\t{}", time_, *p);
    }
}

std::string format_as(const Tissue& x) {
    std::string buffer;
    for (const auto& p: x.extant_cells_) {
        p->format_to_back(buffer);
    }
    return buffer;
}

} // namespace tumopp
