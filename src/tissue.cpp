/*! @file tissue.cpp
    @brief Implementation of Tissue class
*/
#include "tissue.hpp"

#include <wtl/random.hpp>
#include <wtl/iostr.hpp>
#include <wtl/debug.hpp>
#include <wtl/math.hpp>
#include <wtl/numeric.hpp>
#include <wtl/algorithm.hpp>
#include <wtl/genetic.hpp>
#include <wtl/cluster.hpp>
#include <sfmt.hpp>

namespace tumopp {

Tissue::Tissue(
  const size_t initial_size,
  const unsigned dimensions,
  const std::string& coordinate,
  const std::string& local_density_effect,
  const std::string& displacement_path,
  const EventRates& init_event_rates):
  snapshots_(wtl::make_oss()),
  drivers_(wtl::make_oss()) {HERE;
    specimens_.reserve(initial_size * 2u);
    init_coord(dimensions, coordinate);
    init_insert_function(local_density_effect, displacement_path);
    const auto initial_coords = coord_func_->sphere(initial_size);
    const auto origin = std::make_shared<Cell>(
      initial_coords[0], ++id_tail_,
      std::make_shared<EventRates>(init_event_rates)
    );
    tumor_.insert(origin);
    queue_push(origin);
    while (tumor_.size() < initial_size) {
        for (const auto& mother: tumor_) {
            const auto daughter = std::make_shared<Cell>(*mother);
            const auto ancestor = std::make_shared<Cell>(*mother);
            ancestor->set_time_of_death(0.0);
            specimens_.emplace_back(ancestor);
            mother->set_time_of_birth(0.0, ++id_tail_, ancestor);
            daughter->set_time_of_birth(0.0, ++id_tail_, ancestor);
            daughter->set_coord(initial_coords[tumor_.size()]);
            tumor_.insert(daughter);
            // queue_push(mother);
            queue_push(daughter);
            if (tumor_.size() >= initial_size) break;
        }
    }
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
    if (recording_early_growth > 0u) {write_snapshot();}
    bool success = false;
    size_t i = 0;
    double time_snapshot = i_snapshot_ * snapshot_interval;
    while (true) {
        if ((++i % 1000U) == 0U) {DCERR("\r" << tumor_.size());}
        auto it = queue_.begin();
        time_ = it->first;
        if (time_ > max_time || tumor_.size() >= max_size) {
            success = true; // maybe not; but want to exit with record
            break;
        }
        if (time_ > time_snapshot) {
            write_snapshot();
            time_snapshot = ++i_snapshot_ * snapshot_interval;
        }
        const auto mother = std::move(it->second);
        queue_.erase(it);
        if (mother->next_event() == Event::birth) {
            const auto daughter = std::make_shared<Cell>(*mother);
            if (insert(daughter)) {
                const auto ancestor = std::make_shared<Cell>(*mother);
                ancestor->set_time_of_death(time_);
                specimens_.emplace_back(ancestor);
                mother->set_time_of_birth(time_, ++id_tail_, ancestor);
                daughter->set_time_of_birth(time_, ++id_tail_, ancestor);
                drivers_ << mother->mutate();
                drivers_ << daughter->mutate();
                if (tumor_.size() > mutation_timing) {  // once
                    mutation_timing = std::numeric_limits<size_t>::max();
                    drivers_ << daughter->force_mutate();
                }
                queue_push(mother);
                queue_push(daughter);
            } else {
                queue_push(mother);
                continue;  // skip write()
            }
        } else if (mother->next_event() == Event::death) {
            mother->set_time_of_death(time_);
            specimens_.emplace_back(mother);
            tumor_.erase(mother);
            if (tumor_.empty()) break;
        } else {
            migrate(mother);
            queue_push(mother);
        }
        if (tumor_.size() < recording_early_growth) {
            write_snapshot();
        } else {
            recording_early_growth = 0u;  // prevent restart by cell death
        }
    }
    DCERR("\r" << tumor_.size() << std::endl);
    return success;
}

void Tissue::plateau(const double time) {HERE;
    queue_.clear();
    for (auto& p: tumor_) {
        p->increase_death_rate();
        p->set_elapsed(0.0);
        queue_push(p);
    }
    grow(std::numeric_limits<size_t>::max(), time_ + time);
}

void Tissue::treatment(const double death_prob, const size_t num_resistant_cells) {HERE;
    const size_t original_size = tumor_.size();
    std::vector<std::shared_ptr<Cell>> cells;
    cells.reserve(original_size);
    for (const auto& p: queue_) { // for reproducibility
        cells.emplace_back(p.second);
    }
    std::shuffle(cells.begin(), cells.end(), wtl::sfmt64());
    for (size_t i=0; i<original_size; ++i) {
        const auto& p = cells[i];
        if (i >= num_resistant_cells) {
            p->set_cycle_dependent_death(death_prob);
        }
    }
    const size_t margin = 10u * num_resistant_cells + 10u;
    grow(original_size + margin, std::numeric_limits<double>::max());
}

void Tissue::queue_push(const std::shared_ptr<Cell>& x) {
    double dt = x->delta_time(positional_value(x->coord()));
    queue_.emplace_hint(queue_.end(), dt += time_, x);
}

void Tissue::init_insert_function(const std::string& local_density_effect, const std::string& displacement_path) {
    using func_t = std::function<bool(const std::shared_ptr<Cell>&)>;
    using map_sf = std::unordered_map<std::string, func_t>;
    std::unordered_map<std::string, map_sf> swtch;

    swtch["const"].emplace("random", [this](const std::shared_ptr<Cell>& daughter) {
        push(daughter, coord_func_->random_direction(wtl::sfmt64()));
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
        stroll(daughter, coord_func_->random_direction(wtl::sfmt64()));
        return true;
    });
    swtch["step"].emplace("random", [this](const std::shared_ptr<Cell>& daughter) {
        if (num_empty_neighbors(daughter->coord()) == 0U) {return false;}
        push(daughter, coord_func_->random_direction(wtl::sfmt64()));
        return true;
    });
    swtch["step"].emplace("mindrag", [this](const std::shared_ptr<Cell>& daughter) {
        return insert_adjacent(daughter);
    });
    swtch["linear"].emplace("random", [this](const std::shared_ptr<Cell>& daughter) {
        const double prob = proportion_empty_neighbors(daughter->coord());
        if (wtl::generate_canonical(wtl::sfmt64()) < prob) {
            push(daughter, coord_func_->random_direction(wtl::sfmt64()));
            return true;
        }
        return false;
    });
    swtch["linear"].emplace("mindrag", [this](const std::shared_ptr<Cell>& daughter) {
        daughter->set_coord(coord_func_->random_neighbor(daughter->coord(), wtl::sfmt64()));
        return tumor_.insert(daughter).second;
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

void Tissue::push(std::shared_ptr<Cell> moving, const std::valarray<int>& direction) {
    do {
        moving->set_coord(moving->coord() + direction);
    } while (swap_existing(&moving));
}

void Tissue::push_minimum_drag(std::shared_ptr<Cell> moving) {
    do {
        moving->set_coord(moving->coord() + to_nearest_empty(moving->coord()));
    } while (swap_existing(&moving));
}

void Tissue::stroll(std::shared_ptr<Cell> moving, const std::valarray<int>& direction) {
    while (!insert_adjacent(moving)) {
        moving->set_coord(moving->coord() + direction);
        swap_existing(&moving);
    }
}

bool Tissue::insert_adjacent(const std::shared_ptr<Cell>& moving) {
    const auto present_coord = moving->coord();
    auto neighbors = coord_func_->neighbors(present_coord);
    std::shuffle(neighbors.begin(), neighbors.end(), wtl::sfmt64());
    for (auto& x: neighbors) {
        moving->set_coord(x);
        if (tumor_.insert(moving).second) {
            return true;
        }
    }
    moving->set_coord(present_coord);
    return false;
}

bool Tissue::swap_existing(std::shared_ptr<Cell>* x) {
    // The cell must not be in tumor_.
    auto result = tumor_.insert(*x);
    if (result.second) {
        return false;
    } else {
        std::shared_ptr<Cell> existing = std::move(*result.first);
        tumor_.erase(result.first);
        tumor_.insert(std::move(*x));
        x->swap(existing);
        return true;
    }
}

void Tissue::migrate(const std::shared_ptr<Cell>& moving) {
    tumor_.erase(moving);
    const auto orig_pos = moving->coord();
    moving->set_coord(coord_func_->random_neighbor(moving->coord(), wtl::sfmt64()));
    auto result = tumor_.insert(moving);
    if (!result.second) {
        std::shared_ptr<Cell> existing = std::move(*result.first);
        tumor_.erase(result.first);
        tumor_.insert(std::move(moving));
        existing->set_coord(orig_pos);
        tumor_.insert(existing);
    }
}

size_t Tissue::steps_to_empty(std::valarray<int> current, const std::valarray<int>& direction) const {
    size_t steps = 0;
    const auto key = std::make_shared<Cell>();
    do {
        key->set_coord(current += direction);
        ++steps;
    } while (tumor_.find(key) != tumor_.end());
    return steps;
}

std::valarray<int> Tissue::to_nearest_empty(const std::valarray<int>& current, const unsigned search_max) const {
    size_t least_steps = std::numeric_limits<size_t>::max();
    std::valarray<int> best_direction;
    auto directions = coord_func_->directions();
    std::shuffle(directions.begin(), directions.end(), wtl::sfmt64());
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

std::valarray<int> Tissue::roulette_direction(const std::valarray<int>& current) const {
    auto directions = coord_func_->directions();
    std::shuffle(directions.begin(), directions.end(), wtl::sfmt64());
    std::vector<double> roulette;
    for (const auto& d: directions) {
        const auto l = steps_to_empty(current, d);
        if (l == 0U) {return d;}
        roulette.push_back(1.0 / l);
    }
    return directions[wtl::roulette_select(roulette, wtl::sfmt64())];
}

uint_fast8_t Tissue::num_empty_neighbors(const std::valarray<int>& coord) const {
    uint_fast8_t cnt = 0;
    std::shared_ptr<Cell> nb = std::make_shared<Cell>();
    for (const auto& d: coord_func_->directions()) {
        nb->set_coord(coord + d);
        if (tumor_.find(nb) == tumor_.end()) {++cnt;}
    }
    return cnt;
}

std::vector<unsigned> Tissue::generate_neutral_mutations() const {
    std::poisson_distribution<unsigned> poisson(Cell::MUTATION_RATE() * id_tail_);
    const unsigned num_mutants = poisson(wtl::sfmt64());
    std::uniform_int_distribution<unsigned> uniform(1, id_tail_);
    std::vector<unsigned> mutants;
    if (Cell::HAS_AT_LEAST_1_MUTATION_PER_DIVISION()) {
        mutants.reserve(id_tail_ + num_mutants);
        for (unsigned i=1; i<=id_tail_; ++i) mutants.push_back(i);
    } else {
        mutants.reserve(num_mutants);
    }
    for (unsigned i=0; i<num_mutants; ++i) {
        mutants.push_back(uniform(wtl::sfmt64()));
    }
    return mutants;
}

std::ostream& Tissue::write_segsites(std::ostream& ost, const std::vector<std::shared_ptr<Cell>>& samples) const {HERE;
    const size_t sample_size = samples.size();
    const auto mutants = generate_neutral_mutations();
    std::vector<std::vector<unsigned>> flags;
    flags.reserve(sample_size);
    for (const auto& cell: samples) {
        flags.emplace_back(cell->has_mutations_of(mutants));
    }
    flags = wtl::transpose(flags);
    std::vector<std::vector<unsigned>> segsites;
    segsites.reserve(flags.size());
    for (size_t i=0; i<flags.size(); ++i) {
        const auto daf = wtl::sum(flags[i]);
        if (0U < daf && daf < sample_size) segsites.push_back(flags[i]);
    }
    const size_t s = segsites.size();
    ost << "\n//\nsegsites: " << s << "\n";
    if (s > 0U) {
        segsites = wtl::transpose(segsites);
        ost << "positions: ";
        wtl::join(std::vector<int>(s), ost, " ") << "\n";
        for (const auto& x: segsites) {
            wtl::join(x, ost, "") << "\n";
        }
    } else {ost << "\n";}
    return ost;
}

std::vector<std::shared_ptr<Cell>>
Tissue::sample_bulk(const std::shared_ptr<Cell>& center, const size_t n) const {
    std::multimap<double, std::shared_ptr<Cell>> ordered;
    for (const auto& p: tumor_) {
        ordered.emplace(coord_func_->euclidean_distance(p->coord() - center->coord()), p);
    }
    std::vector<std::shared_ptr<Cell>> sampled;
    sampled.reserve(n);
    for (const auto& p: ordered) {
        sampled.emplace_back(p.second);
        if (sampled.size() >= n) break;
    }
    return sampled;
}

std::vector<std::shared_ptr<Cell>>
Tissue::sample_medoids(const size_t n) const {HERE;
    std::vector<std::shared_ptr<Cell>> cells(tumor_.begin(), tumor_.end());
    std::vector<std::valarray<double>> points;
    points.reserve(n);
    for (const auto& p: cells) {
        points.emplace_back(coord_func_->continuous(p->coord()));
    }
    auto clusters = wtl::cluster::pam(points, n, wtl::sfmt64());
    std::vector<std::shared_ptr<Cell>> sampled;
    sampled.reserve(n);
    for (const auto i: clusters.medoids()) {
        sampled.emplace_back(cells[i]);
    }
    return sampled;
}

std::vector<std::shared_ptr<Cell>> Tissue::sample_random(const size_t n) const {HERE;
    return wtl::sample(std::vector<std::shared_ptr<Cell>>(tumor_.begin(), tumor_.end()), n, wtl::sfmt64());
}

std::vector<std::shared_ptr<Cell>> Tissue::sample_section(const size_t n) const {HERE;
    std::vector<std::shared_ptr<Cell>> section;
    section.reserve(static_cast<size_t>(coord_func_->cross_section(tumor_.size())));
    for (const auto& p: tumor_) {
        if (p->coord()[2] == 0) {section.push_back(p);}
    }
    return wtl::sample(section, n, wtl::sfmt64());
}

std::string Tissue::pairwise_distance(const size_t npair) const {HERE;
    auto oss = wtl::make_oss(6);
    oss << "genealogy\tgraph\teuclidean\n";
    if (npair == 0u) return oss.str();
    auto samples = sample_random(2 * npair);
    std::shuffle(samples.begin(), samples.end(), wtl::sfmt64());
    //TODO: should be randam sampling from all possible pairs
    const auto end = samples.cend();
    for (auto it=samples.cbegin(); it!=end; ++it) {
        const auto& lhs = *(*it);
        const auto& rhs = *(*(++it));
        const auto diff = lhs.coord() - rhs.coord();
        oss << lhs.branch_length(rhs) << "\t"
            << coord_func_->graph_distance(diff) << "\t"
            << coord_func_->euclidean_distance(diff) << "\n";
    }
    return oss.str();
}

std::string Tissue::header() const {
    return Cell::header() + "\tphi\n";
}

void Tissue::write(std::ostream& ost, const Cell& cell) const {
    cell.write(ost) << "\t"
       << static_cast<unsigned>(num_empty_neighbors(cell.coord())) << "\n";
}

std::string Tissue::specimens() {
    for (const auto& p: tumor_) {
        specimens_.emplace_back(p);
    }
    tumor_.clear();
    queue_.clear();
    std::ostringstream oss;
    oss << header();
    for (const auto& p: specimens_) {
        write(oss, *p);
    }
    return oss.str();
}

void Tissue::write_snapshot() {
    for (const auto& p: tumor_) {
        write(snapshots_ << time_ << "\t", *p);
    }
}

//! Stream operator for debug print
std::ostream& operator<< (std::ostream& ost, const Tissue& tissue) {
    for (const auto& p: tissue.tumor_) {
        ost << *p << "\n";
    }
    return ost;
}

} // namespace tumopp
