// -*- mode: c++; coding: utf-8 -*-
/*! @file tissue.cpp
    @brief Implementation of Tissue class
*/
#include "tissue.hpp"

#include <cstdlib>
#include <cassert>
#include <iostream>
#include <sstream>

#include <boost/program_options.hpp>

#include <cxxwtils/prandom.hpp>
#include <cxxwtils/iostr.hpp>
#include <cxxwtils/gz.hpp>
#include <cxxwtils/debug.hpp>
#include <cxxwtils/algorithm.hpp>


size_t Tissue::DIMENSIONS_ = 3;
std::string Tissue::COORDINATE_ = "moore";
std::string Tissue::PACKING_ = "push";
double Tissue::SIGMA_E_ = std::numeric_limits<double>::infinity();
size_t Tissue::INITIAL_SIZE_ = 1;
size_t Tissue::MAX_SIZE_ = 16384;

//! Program options
/*! @return Program options description

    Command line option | Symbol         | Variable
    --------------------| -------------- | -------------------------
    `-D,--dimensions`   | -              | Tissue::DIMENSIONS_
    `-C,--coord`        | -              | Tissue::COORDINATE_
    `-P,--packing`      | -              | Tissue::PACKING_
    `-g,--peripheral`   | \f$\sigma_E\f$ | Tissue::SIGMA_E_
    `-O,--origin`       | -              | Tissue::INITIAL_SIZE_
    `-N,--max`          | -              | Simulation::MAX_SIZE
*/
boost::program_options::options_description& Tissue::opt_description() {
    namespace po = boost::program_options;
    static po::options_description desc{"Tissue"};
    desc.add_options()
        ("dimensions,D", po::value<size_t>(&DIMENSIONS_)->default_value(DIMENSIONS_))
        ("coord,C", po::value<std::string>(&COORDINATE_)->default_value(COORDINATE_))
        ("packing,P", po::value<std::string>(&PACKING_)->default_value(PACKING_))
        ("peripheral,g", po::value<double>(&SIGMA_E_)->default_value(SIGMA_E_))
        ("origin,O", po::value<size_t>(&INITIAL_SIZE_)->default_value(INITIAL_SIZE_))
        ("max,N", po::value<size_t>(&MAX_SIZE_)->default_value(MAX_SIZE_))
    ;
    return desc;
}

Tissue::Tissue() {
    if (COORDINATE_ == "neumann") {coord_func_ = std::make_unique<Neumann>(DIMENSIONS_);}
    else if (COORDINATE_ == "moore") {coord_func_ = std::make_unique<Moore>(DIMENSIONS_);}
    else if (COORDINATE_ == "hex") {coord_func_ = std::make_unique<Hexagonal>(DIMENSIONS_);}
    specimens_.precision(std::numeric_limits<double>::max_digits10);
    snapshots_.precision(std::numeric_limits<double>::max_digits10);
    for (const auto& coord: coord_func_->sphere(INITIAL_SIZE_)) {
        auto x = std::make_shared<Cell>(coord, ++id_tail_);
        tumor_.insert(x);
        queue_push(x->delta_time(positional_value(x->coord())), x);
    }
}

bool Tissue::grow(const size_t max_size) {HERE;
    snap(snapshots_);
    bool taking_snapshots = true;
    while (tumor_.size() < max_size) {
        auto it = queue_.begin();
        time_ = it->first;
        const auto mother = it->second;
        queue_.erase(it);
        if (mother->is_dividing()) {
            const auto daughter = std::make_shared<Cell>(*mother);
            if (insert(daughter)) {
                mother->set_time_of_death(time_);
                collect(specimens_, *mother);
                mother->set_time_of_death(0.0);
                mother->set_time_of_birth(time_, ++id_tail_);
                daughter->set_time_of_birth(time_, ++id_tail_);
                mother->mutate();
                daughter->mutate();
                queue_push(mother->delta_time(positional_value(mother->coord())), mother);
                queue_push(daughter->delta_time(positional_value(daughter->coord())), daughter);
            } else {
                queue_push(mother->delta_time(positional_value(mother->coord())), mother);
                continue;  // skip snap()
            }
        } else if (mother->is_dying()) {
            mother->set_time_of_death(time_);
            collect(specimens_, *mother);
            tumor_.erase(mother);
            if (tumor_.empty()) return false;
        } else {
            migrate(mother);
            queue_push(mother->delta_time(positional_value(mother->coord())), mother);
        }
        if (taking_snapshots && tumor_.size() < 128) {
            snap(snapshots_);
        } else {
            taking_snapshots = false;  // prevent restart by cell death
        }
    }
    snap(specimens_);
    return true;
}

void Tissue::queue_push(double delta_t, const std::shared_ptr<Cell>& x) {
    queue_.insert(queue_.end(), std::make_pair(delta_t += time_, x));
}

bool Tissue::insert(const std::shared_ptr<Cell>& daughter) {
    if (PACKING_ == "push") {
        push(daughter, coord_func_->random_direction(wtl::sfmt()));
    } else if (PACKING_ == "pushn") {
        push(daughter, to_nearest_empty(daughter->coord()));
    } else if (PACKING_ == "pushne") {
        pushn_everytime(daughter);
    } else if (PACKING_ == "fillpush") {
        fill_push(daughter, coord_func_->random_direction(wtl::sfmt()));
    } else if (PACKING_ == "fill") {  //! @todo incorrect time scale
        return fill_empty(daughter);
    } else if (PACKING_ == "empty") {
        return insert_neighbor(daughter);
    } else {
        std::cout << FILE_LINE_PRETTY
            << "\nERROR: unknown option value --packing=" << PACKING_
            << std::endl;
        std::abort();
    }
    return true;
}

void Tissue::push(std::shared_ptr<Cell> moving, const std::vector<int>& direction) {
    do {
        moving->set_coord(moving->coord() + direction);
    } while (swap_existing(&moving));
}

void Tissue::pushn_everytime(std::shared_ptr<Cell> moving) {
    do {
        moving->set_coord(moving->coord() + to_nearest_empty(moving->coord()));
    } while (swap_existing(&moving));
}

void Tissue::fill_push(std::shared_ptr<Cell> moving, const std::vector<int>& direction) {
    while (!fill_empty(moving)) {
        moving->set_coord(moving->coord() + direction);
        swap_existing(&moving);
    }
}

bool Tissue::fill_empty(const std::shared_ptr<Cell>& moving) {
    const auto present_coord = moving->coord();
    auto neighbors = coord_func_->neighbors(present_coord);
    std::shuffle(neighbors.begin(), neighbors.end(), wtl::sfmt());
    for (auto& x: neighbors) {
        moving->set_coord(x);
        if (tumor_.insert(moving).second) {
            // end if found and filled an empty space
            return true;
        }
    }
    moving->set_coord(present_coord);
    return false;
}

bool Tissue::insert_neighbor(const std::shared_ptr<Cell>& daughter) {
    daughter->set_coord(coord_func_->random_neighbor(daughter->coord(), wtl::sfmt()));
    return tumor_.insert(daughter).second;
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

size_t Tissue::steps_to_empty(std::vector<int> current, const std::vector<int>& direction) const {
    size_t steps = 0;
    const auto key = std::make_shared<Cell>();
    do {
        key->set_coord(current += direction);
        ++steps;
    } while (tumor_.find(key) != tumor_.end());
    return steps;
}

std::vector<int> Tissue::to_nearest_empty(const std::vector<int>& current, size_t search_max) const {
    search_max = std::min(search_max, coord_func_->directions().size());
    size_t least_steps = std::numeric_limits<size_t>::max();
    std::vector<int> best_direction;
    for (const auto& d: wtl::sample_knuth(coord_func_->directions(), search_max, wtl::sfmt())) {
        auto n = steps_to_empty(current, d);
        if (n < least_steps) {
            least_steps = n;
            best_direction = d;
        }
    }
    return best_direction;
}

std::vector<std::vector<int>> Tissue::empty_neighbors(const std::vector<int>& coord) const {
    std::vector<std::vector<int>> output;
    std::shared_ptr<Cell> nb = std::make_shared<Cell>();
    for (const auto& d: coord_func_->directions()) {
        nb->set_coord(coord + d);
        if (tumor_.find(nb) != tumor_.end()) {
            output.push_back(nb->coord());
        }
    }
    return output;
}

double Tissue::positional_value(const std::vector<int>& coord) const {
    if (SIGMA_E_ > 1e9 | tumor_.size() <= 8) return 1.0;
    double rel_d = coord_func_->euclidean_distance(coord)
                   / coord_func_->radius(tumor_.size());
    double exponent = wtl::pow(std::max(0.0, 1.0 - rel_d), 2);
    exponent /= SIGMA_E_;
    return std::exp(-exponent);
}

std::vector<size_t> Tissue::generate_neutral_mutations() const {
    std::poisson_distribution<size_t> poisson(Cell::MUTATION_RATE() * id_tail_);
    const size_t num_mutants = poisson(wtl::sfmt());
    std::uniform_int_distribution<size_t> uniform(0, id_tail_);
    std::vector<size_t> mutants;
    mutants.reserve(num_mutants);
    for (size_t i=0; i<num_mutants; ++i) {
        mutants.push_back(uniform(wtl::sfmt()));
    }
    return mutants;
}

std::ostream& Tissue::write_segsites(std::ostream& ost, const std::vector<std::shared_ptr<Cell>>& samples) const {HERE;
    const auto mutants = generate_neutral_mutations();
    std::vector<std::vector<int>> flags;
    flags.reserve(samples.size());
    for (const auto& cell: samples) {
        flags.push_back(cell->has_mutations_of(mutants));
    }
    wtl::transpose(&flags);
    std::vector<std::vector<int>> segsites;
    segsites.reserve(flags.size());
    for (size_t i=0; i<flags.size(); ++i) {
        if (wtl::sum(flags[i]) > 0) segsites.push_back(flags[i]);
    }
    const size_t s = segsites.size();
    wtl::transpose(&segsites);
    ost << "\n//\nsegsites: " << s << "\n";
    if (s > 0) {
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
    section.reserve(coord_func_->cross_section(MAX_SIZE_));
    for (const auto& p: tumor_) {
        if (p->coord()[2] == 0) {section.push_back(p);}
    }
    return wtl::sample(section, n, wtl::sfmt());
}

std::string Tissue::header() const {HERE;
    std::ostringstream oss;
    oss.precision(std::numeric_limits<double>::max_digits10);
    oss << "time" << sep_ << Cell::header(DIMENSIONS_, sep_);
    return oss.str();
}

void Tissue::collect(std::ostream& ost, const Cell& cell) {
    cell.write(ost << time_ << sep_, sep_);
}

void Tissue::snap(std::ostream& ost) {
    for (const auto& p: tumor_) {
        collect(ost, *p);
    }
}

//! Stream operator for debug print
std::ostream& operator<< (std::ostream& ost, const Tissue& tissue) {
    for (const auto& p: tissue.tumor_) {
        ost << *p << "\n";
    }
    return ost;
}

template <class T> inline
void test_coordinate(const std::vector<int>& v) {HERE;
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

template <class T> inline
void test_radius() {HERE;
    Tissue tissue;
    tissue.set_coord<T>();
    tissue.grow(100);
    std::cerr << tissue.size() << ": " << tissue.radius() << std::endl;
    tissue.grow(1000);
    std::cerr << tissue.size() << ": " << tissue.radius() << std::endl;
}

void Tissue::unit_test() {HERE;
    std::cerr.precision(15);

    Tissue tissue;
    tissue.grow(10);
    std::cerr << tissue << std::endl;
    std::cerr << tissue.header();
    tissue.snap(std::cerr);

    const std::vector<int> v2{3, -2};
    test_coordinate<Neumann>(v2);
    test_coordinate<Moore>(v2);
    test_coordinate<Hexagonal>(v2);

    const std::vector<int> v3{3, -2, 1};
    test_coordinate<Neumann>(v3);
    test_coordinate<Moore>(v3);
    test_coordinate<Hexagonal>(v3);

    test_radius<Moore>();
    test_radius<Hexagonal>();
}
