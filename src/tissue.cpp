// -*- mode: c++; coding: utf-8 -*-
/*! @file tissue.cpp
    @brief Implementation of Tissue class
*/
#include "tissue.hpp"

#include <sfmt.hpp>
#include <wtl/prandom.hpp>
#include <wtl/iostr.hpp>
#include <wtl/debug.hpp>
#include <wtl/math.hpp>
#include <wtl/numeric.hpp>
#include <wtl/algorithm.hpp>
#include <wtl/genetic.hpp>

#include <boost/program_options.hpp>

#include <iostream>
#include <sstream>

namespace tumopp {

unsigned int Tissue::DIMENSIONS_;
std::string Tissue::COORDINATE_;
std::string Tissue::LOCAL_DENSITY_EFFECT_;
std::string Tissue::DISPLACEMENT_PATH_;
double Tissue::SIGMA_E_;
size_t Tissue::INITIAL_SIZE_;
size_t Tissue::RECORDING_EARLY_GROWTH_;
size_t Tissue::MUTATION_TIMING_;

//! Program options
/*! @ingroup params
    @return Program options description

    Command line option | Symbol         | Variable
    ------------------- | -------------- | -------------------------------
    `-D,--dimensions`   | -              | Tissue::DIMENSIONS_
    `-C,--coord`        | -              | Tissue::COORDINATE_
    `-L,--local`        | \f$E_2\f$      | Tissue::LOCAL_DENSITY_EFFECT_
    `-P,--path`         | -              | Tissue::DISPLACEMENT_PATH_
    `-g,--peripheral`   | \f$\sigma_E\f$ | Tissue::SIGMA_E_
    `-O,--origin`       | \f$N_0\f$      | Tissue::INITIAL_SIZE_
    `-R,--record`       | -              | Tissue::RECORDING_EARLY_GROWTH_
    `-U,--mutate`       | \f$N_\mu\f$    | Tissue::MUTATION_TIMING_
*/
boost::program_options::options_description Tissue::opt_description() {
    namespace po = boost::program_options;
    po::options_description desc{"Tissue"};
    desc.add_options()
        ("dimensions,D", po::value(&DIMENSIONS_)->default_value(3))
        ("coord,C", po::value(&COORDINATE_)->default_value("moore"))
        ("local,L", po::value(&LOCAL_DENSITY_EFFECT_)->default_value("const"))
        ("path,P", po::value(&DISPLACEMENT_PATH_)->default_value("default"))
        ("peripheral,g", po::value(&SIGMA_E_)->default_value(std::numeric_limits<double>::infinity()))
        ("origin,O", po::value(&INITIAL_SIZE_)->default_value(1))
        ("record,R", po::value(&RECORDING_EARLY_GROWTH_)->default_value(0))
        ("mutate,U", po::value(&MUTATION_TIMING_)->default_value(std::numeric_limits<size_t>::max()))
    ;
    return desc;
}

void Tissue::init() {HERE;
    tumor_.clear();
    queue_.clear();
    time_ = 0.0;
    id_tail_ = 0;
    (specimens_ = wtl::make_oss()) << header();
    (snapshots_ = wtl::make_oss()) << header();
    (drivers_ = wtl::make_oss()) << "id\ttype\tcoef\n";
    init_coord();
    init_insert_function();
    const auto initial_coords = coord_func_->sphere(INITIAL_SIZE_);
    const auto origin = std::make_shared<Cell>(initial_coords[0], ++id_tail_);
    tumor_.insert(origin);
    queue_push(origin);
    while (tumor_.size() < INITIAL_SIZE_) {
        for (const auto& mother: tumor_) {
            const auto daughter = std::make_shared<Cell>(*mother);
            mother->set_time_of_death(0.0);
            write(specimens_, *mother);
            mother->set_time_of_birth(0.0, ++id_tail_);
            daughter->set_time_of_birth(0.0, ++id_tail_);
            daughter->set_coord(initial_coords[tumor_.size()]);
            tumor_.insert(daughter);
            queue_push(daughter);
            if (tumor_.size() >= INITIAL_SIZE_) break;
        }
    }
    write(snapshots_);
}

void Tissue::init_coord() {HERE;
    std::unordered_map<std::string, std::unique_ptr<Coord>> swtch;
    swtch["neumann"] = std::make_unique<Neumann>(DIMENSIONS_);
    swtch["moore"] = std::make_unique<Moore>(DIMENSIONS_);
    swtch["hex"] = std::make_unique<Hexagonal>(DIMENSIONS_);
    try {
        coord_func_ = std::move(swtch.at(COORDINATE_));
    } catch (std::exception& e) {
        std::ostringstream oss;
        oss << std::endl << FILE_LINE_PRETTY
            << "\nInvalid value for -C (" << COORDINATE_ << "); choose from "
            << wtl::keys(swtch);
        throw std::runtime_error(oss.str());
    }
}

bool Tissue::grow(const size_t max_size, const double plateau_time) {HERE;
    init();
    bool success = false;
    double time_stopped = 0.0;
    const size_t max_time = static_cast<size_t>(std::log2(max_size) * 10);
    size_t i = 0;
    while (true) {
        if (tumor_.size() >= max_size && time_stopped == 0.0) {
            // only the first time that tumor_.size() has reached max_size
            success = true;
            time_stopped = time_;
            queue_.clear();
            for (auto& p: tumor_) {
                p->increase_death_rate();
                queue_push(p);
            }
        }
        if ((++i % 256U) == 0U) {DCERR("\r" << tumor_.size());}
        auto it = queue_.begin();
        time_ = it->first;
        if (time_stopped > 0.0 && time_ - time_stopped > plateau_time) {break;}
        if (time_ > max_time) {success = true; break;}
        const auto mother = it->second;
        queue_.erase(it);
        if (mother->next_event() == Event::birth) {
            const auto daughter = std::make_shared<Cell>(*mother);
            if (insert(daughter)) {
                mother->set_time_of_death(time_);
                write(specimens_, *mother);
                mother->set_time_of_death(0.0);
                mother->set_time_of_birth(time_, ++id_tail_);
                daughter->set_time_of_birth(time_, ++id_tail_);
                drivers_ << mother->mutate();
                drivers_ << daughter->mutate();
                if (tumor_.size() > MUTATION_TIMING_) {  // once
                    MUTATION_TIMING_ = std::numeric_limits<size_t>::max();
                    drivers_ << daughter->force_mutate();
                }
                queue_push(mother);
                queue_push(daughter);
            } else {
                if (mother->frustration() > 250U) break;
                queue_push(mother);
                continue;  // skip write()
            }
        } else if (mother->next_event() == Event::death) {
            mother->set_time_of_death(time_);
            write(specimens_, *mother);
            tumor_.erase(mother);
            if (tumor_.empty()) break;
        } else {
            migrate(mother);
            queue_push(mother);
        }
        if (tumor_.size() < RECORDING_EARLY_GROWTH_) {
            write(snapshots_);
        } else {
            RECORDING_EARLY_GROWTH_ = 0;  // prevent restart by cell death
        }
    }
    DCERR(std::endl);
    write(specimens_);
    return success;
}

void Tissue::queue_push(const std::shared_ptr<Cell>& x) {
    double dt = x->delta_time(positional_value(x->coord()));
    queue_.emplace_hint(queue_.end(), dt += time_, x);
}

void Tissue::init_insert_function() {
    typedef std::function<bool(const std::shared_ptr<Cell>&)> func_t;
    typedef std::unordered_map<std::string, func_t> map_sf;
    std::unordered_map<std::string, map_sf> swtch;

    swtch["const"]["random"] = [this](const std::shared_ptr<Cell>& daughter) {
        push(daughter, coord_func_->random_direction(wtl::sfmt()));
        return true;
    };
    swtch["const"]["mindrag"] = [this](const std::shared_ptr<Cell>& daughter) {
        push_minimum_drag(daughter);
        return true;
    };
    swtch["const"]["minstraight"] = [this](const std::shared_ptr<Cell>& daughter) {
        push(daughter, to_nearest_empty(daughter->coord()));
        return true;
    };
    swtch["const"]["roulette"] = [this](const std::shared_ptr<Cell>& daughter) {
        push(daughter, roulette_direction(daughter->coord()));
        return true;
    };
    swtch["const"]["stroll"] = [this](const std::shared_ptr<Cell>& daughter) {
        stroll(daughter, coord_func_->random_direction(wtl::sfmt()));
        return true;
    };
    swtch["step"]["random"] = [this](const std::shared_ptr<Cell>& daughter) {
        if (num_empty_neighbors(daughter->coord()) == 0U) {return false;}
        push(daughter, coord_func_->random_direction(wtl::sfmt()));
        return true;
    };
    swtch["step"]["mindrag"] = [this](const std::shared_ptr<Cell>& daughter) {
        return insert_adjacent(daughter);
    };
    swtch["linear"]["random"] = [this](const std::shared_ptr<Cell>& daughter) {
        const double prob = proportion_empty_neighbors(daughter->coord());
        if (!std::bernoulli_distribution(prob)(wtl::sfmt())) {return false;}
        push(daughter, coord_func_->random_direction(wtl::sfmt()));
        return true;
    };
    swtch["linear"]["mindrag"] = [this](const std::shared_ptr<Cell>& daughter) {
        daughter->set_coord(coord_func_->random_neighbor(daughter->coord(), wtl::sfmt()));
        return tumor_.insert(daughter).second;
    };
    swtch["const"]["default"] = swtch["const"]["random"];
    swtch["step"]["default"] = swtch["step"]["mindrag"];
    swtch["linear"]["default"] = swtch["linear"]["mindrag"];
    try {
        insert = swtch.at(LOCAL_DENSITY_EFFECT_).at(DISPLACEMENT_PATH_);
    } catch (std::exception& e) {
        std::ostringstream oss;
        oss << std::endl << FILE_LINE_PRETTY
            << "\nInvalid value for -L/-P ("
            << LOCAL_DENSITY_EFFECT_ << "/" << DISPLACEMENT_PATH_
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
    std::shuffle(neighbors.begin(), neighbors.end(), wtl::sfmt());
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
    moving->set_coord(coord_func_->random_neighbor(moving->coord(), wtl::sfmt()));
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

std::valarray<int> Tissue::to_nearest_empty(const std::valarray<int>& current, const unsigned int search_max) const {
    size_t least_steps = std::numeric_limits<size_t>::max();
    std::valarray<int> best_direction;
    auto directions = coord_func_->directions();
    std::shuffle(directions.begin(), directions.end(), wtl::sfmt());
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
    std::shuffle(directions.begin(), directions.end(), wtl::sfmt());
    std::vector<double> roulette;
    for (const auto& d: directions) {
        const auto l = steps_to_empty(current, d);
        if (l == 0U) {return d;}
        roulette.push_back(1.0 / l);
    }
    return directions[wtl::roulette_select(roulette, wtl::sfmt())];
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

double Tissue::positional_value(const std::valarray<int>& coord) const {
    if ((SIGMA_E_ > 1e9) | (tumor_.size() <= 8U)) return 1.0;
    double rel_d = coord_func_->euclidean_distance(coord)
                   / coord_func_->radius(tumor_.size());
    double exponent = wtl::pow(std::max(0.0, 1.0 - rel_d), 2);
    exponent /= SIGMA_E_;
    return std::exp(-exponent);
}

std::vector<size_t> Tissue::generate_neutral_mutations() const {
    std::poisson_distribution<unsigned int> poisson(Cell::MUTATION_RATE() * id_tail_);
    const unsigned int num_mutants = poisson(wtl::sfmt());
    std::uniform_int_distribution<size_t> uniform(0, id_tail_);
    std::vector<size_t> mutants;
    mutants.reserve(num_mutants);
    for (unsigned int i=0; i<num_mutants; ++i) {
        mutants.push_back(uniform(wtl::sfmt()));
    }
    return mutants;
}

std::ostream& Tissue::write_segsites(std::ostream& ost, const std::vector<std::shared_ptr<Cell>>& samples) const {HERE;
    const size_t sample_size = samples.size();
    const auto mutants = generate_neutral_mutations();
    std::vector<std::vector<unsigned int>> flags;
    flags.reserve(sample_size);
    for (const auto& cell: samples) {
        flags.push_back(cell->has_mutations_of(mutants));
    }
    flags = wtl::transpose(flags);
    std::vector<std::vector<unsigned int>> segsites;
    segsites.reserve(flags.size());
    for (size_t i=0; i<flags.size(); ++i) {
        const auto daf = wtl::sum(flags[i]);
        if (0U < daf && daf < sample_size) segsites.push_back(flags[i]);
    }
    const size_t s = segsites.size();
    ost << "\n//\nsegsites: " << s << "\n";
    if (s > 0U) {
        segsites = wtl::transpose(segsites);
        ost << "positions: "
            << wtl::join(std::vector<int>(s), " ") << "\n";
        for (const auto& x: segsites) {
            ost << wtl::join(x, "") << "\n";
        }
    } else {ost << "\n";}
    return ost;
}

std::vector<std::shared_ptr<Cell>> Tissue::sample_random(const size_t n) const {HERE;
    return wtl::sample(std::vector<std::shared_ptr<Cell>>(tumor_.begin(), tumor_.end()), n, wtl::sfmt());
}

std::vector<std::shared_ptr<Cell>> Tissue::sample_section(const size_t n) const {HERE;
    std::vector<std::shared_ptr<Cell>> section;
    section.reserve(coord_func_->cross_section(tumor_.size()));
    for (const auto& p: tumor_) {
        if (p->coord()[2] == 0) {section.push_back(p);}
    }
    return wtl::sample(section, n, wtl::sfmt());
}

std::string Tissue::pairwise_distance(const size_t npair) const {HERE;
    auto oss = wtl::make_oss(6);
    oss << "genealogy\tgraph\teuclidean\n";
    auto samples = sample_random(2 * npair);
    std::shuffle(samples.begin(), samples.end(), wtl::sfmt());
    //TODO: should be randam sampling from all possible pairs
    const auto end = samples.cend();
    for (auto it=samples.cbegin(); it!=end; ++it) {
        const auto& lhs = *(*it);
        const auto& rhs = *(*(++it));
        const auto diff = lhs.coord() - rhs.coord();
        oss << lhs.branch_length(rhs) << "\t"
            << coord_func()->graph_distance(diff) << "\t"
            << coord_func()->euclidean_distance(diff) << "\n";
    }
    return oss.str();
}

std::string Tissue::header() const {HERE;
    auto oss = wtl::make_oss();
    oss << "time\t" << Cell::header() << "\tphi\n";
    return oss.str();
}

void Tissue::write(std::ostream& ost, const Cell& cell) const {
    cell.write(ost << time_ << "\t") << "\t"
       << static_cast<unsigned int>(num_empty_neighbors(cell.coord())) << "\n";
}

void Tissue::write(std::ostream& ost) const {
    for (const auto& p: tumor_) {
        write(ost, *p);
    }
}

//! Stream operator for debug print
std::ostream& operator<< (std::ostream& ost, const Tissue& tissue) {
    for (const auto& p: tissue.tumor_) {
        ost << *p << "\n";
    }
    return ost;
}

//! Test coordinate system
template <class T> inline
void test_coordinate(const std::valarray<int>& v) {HERE;
    T coord(v.size());
    std::cerr << coord.core() << std::endl;
    std::cerr << coord.sphere(20) << std::endl;
    for (auto x: coord.directions()) {
        std::cerr << x  << ": " << coord.euclidean_distance(x) << std::endl;
    }
    std::cerr << coord.neighbors(v) << std::endl;
    std::cerr << coord.graph_distance(v) << std::endl;
    std::cerr << coord.euclidean_distance(v) << std::endl;
}

//! Test radius function
template <class T> inline
void test_radius() {HERE;
    Tissue tissue;
    tissue.init_coord_test<T>();
    tissue.grow(100);
    std::cerr << tissue.size() << ": " << tissue.radius() << std::endl;
    tissue.grow(1000);
    std::cerr << tissue.size() << ": " << tissue.radius() << std::endl;
}

void Tissue::test() {HERE;
    std::cerr.precision(15);

    Tissue tissue;
    tissue.grow(10);
    std::cerr << tissue << std::endl;
    std::cerr << tissue.header();
    tissue.write(std::cerr);

    const std::valarray<int> v2{3, -2};
    test_coordinate<Neumann>(v2);
    test_coordinate<Moore>(v2);
    test_coordinate<Hexagonal>(v2);

    const std::valarray<int> v3{3, -2, 1};
    test_coordinate<Neumann>(v3);
    test_coordinate<Moore>(v3);
    test_coordinate<Hexagonal>(v3);

    test_radius<Moore>();
    test_radius<Hexagonal>();
}

} // namespace tumopp
