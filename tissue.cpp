// -*- mode: c++; coding: utf-8 -*-
/*! @file tissue.cpp
    @brief Implementation of Tissue class
*/
#include "tissue.h"

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
std::string Tissue::SCHEDULE_ = "random";
std::string Tissue::PACKING_ = "push";

//! Program options
/*! @return Program options description

    Command line option | Symbol | Variable
    --------------------| ------ | -------------------------
    `-D,--dimensions`   | -      | Tissue::DIMENSIONS_
    `-C,--coord`        | -      | Tissue::COORDINATE_
    `-S,--schedule`     | -      | Tissue::SCHEDULE_
    `-P,--packing`      | -      | Tissue::PACKING_
*/
boost::program_options::options_description& Tissue::opt_description() {
    namespace po = boost::program_options;
    static po::options_description desc{"Tissue"};
    desc.add_options()
        ("dimensions,D", po::value<size_t>(&DIMENSIONS_)->default_value(DIMENSIONS_))
        ("coord,C", po::value<std::string>(&COORDINATE_)->default_value(COORDINATE_))
        ("schedule,S", po::value<std::string>(&SCHEDULE_)->default_value(SCHEDULE_))
        ("packing,P", po::value<std::string>(&PACKING_)->default_value(PACKING_))
    ;
    return desc;
}

void Tissue::queue_push(double t, const std::shared_ptr<Cell>& x) {
    queue_.insert(queue_.end(), std::make_pair(t, x));
}

void Tissue::grow(const size_t max_size) {HERE;
    stock_.reserve(max_size);
    if (tumor_.empty()) {
        for (const auto& coord: coord_func_->core()) {
            auto x = std::make_shared<Cell>(coord);
            stock_.push_back(x);
            tumor_.insert(x);
            queue_push(x->delta_time(), x);
        }
    }
    evolution_history_.reserve(max_size);
    evolution_history_.push_back(snapshot());
    while (tumor_.size() < max_size) {
        auto it = queue_.begin();
        const auto mother = it->second;
        if (wtl::prandom().bernoulli(it->second->birth_given_event())) {
            const auto daughter = std::make_shared<Cell>(*mother);
            if (insert(daughter)) {
                daughter->set_time_of_birth(it->first);
                stock_.push_back(daughter);
                if (wtl::prandom().bernoulli(daughter->mutation_rate())) {
                    daughter->mutate();
                    mutation_coords_.push_back(daughter->coord());
                    mutation_stages_.push_back(tumor_.size());
                }
                if (wtl::prandom().bernoulli(mother->mutation_rate())) {
                    mother->mutate();
                    mutation_coords_.push_back(mother->coord());
                    mutation_stages_.push_back(tumor_.size());
                }
                queue_push(it->first + mother->delta_time(), mother);
                queue_push(it->first + daughter->delta_time(), daughter);
            } else {
                queue_push(it->first + mother->delta_time(), mother);
            }
        } else {
            it->second->set_time_of_death(it->first);
            tumor_.erase(it->second);
        }
        queue_.erase(it);
    }
}


void Tissue::grow_old(const size_t max_size) {HERE;
    stock_.reserve(max_size);
    // unit is generation time under BIRTH_RATE_ == 1.0
    double time = 0;
    // sum of birth and death rate, which determines waiting time
    double event_rate = 0;
    if (tumor_.empty()) {
        for (const auto& coord: coord_func_->core()) {
            auto x = std::make_shared<Cell>(coord);
            stock_.push_back(x);
            tumor_.insert(x);
            event_rate += x->instantaneous_event_rate();
        }
    }
    evolution_history_.reserve(max_size);
    evolution_history_.push_back(snapshot());
    while (tumor_.size() < max_size) {
        std::vector<std::shared_ptr<Cell>> mothers;
        if (SCHEDULE_ == "random") {
            auto it = std::next(tumor_.begin(), wtl::prandom().randrange(tumor_.size()));
            mothers.push_back(*it);
        } else if (SCHEDULE_ == "even") {
            mothers.reserve(tumor_.size());
            for (auto it=tumor_.begin(); it!=tumor_.end(); ++it) {
                mothers.push_back(*it);
            }
        } else {exit(1);}
        for (const auto& mother: mothers) {
            const bool dying = wtl::prandom().bernoulli(mother->death_rate());
            bool dividing = wtl::prandom().bernoulli(mother->birth_rate());
            if (dividing) {
                time += 1.0 / event_rate;
                const auto daughter = std::make_shared<Cell>(*mother);  // Cell copy ctor
                dividing = insert(daughter);
                if (dividing) {
                    daughter->set_time_of_birth(time);
                    stock_.push_back(daughter);
                    if (wtl::prandom().bernoulli(daughter->mutation_rate())) {
                        daughter->mutate();
                        mutation_coords_.push_back(daughter->coord());
                        mutation_stages_.push_back(tumor_.size());
                    }
                    if (!dying && wtl::prandom().bernoulli(mother->mutation_rate())) {
                        event_rate -= mother->instantaneous_event_rate();
                        mother->mutate();
                        mutation_coords_.push_back(mother->coord());
                        mutation_stages_.push_back(tumor_.size());
                        event_rate += mother->instantaneous_event_rate();
                    }
                }
            }
            if (dying) {
                time += 1.0 / event_rate;
                mother->set_time_of_death(time);
                tumor_.erase(mother);
                event_rate -= mother->instantaneous_event_rate();
            }
            // event_rate should be changed after time increment;
            if (dividing) {
                event_rate += stock_.back()->instantaneous_event_rate();
            }
            if (tumor_.size() <= 128 && (dividing || dying)) {
                evolution_history_.push_back(snapshot());
            }
        }
    }
}

bool Tissue::insert(const std::shared_ptr<Cell>& daughter) {
    if (PACKING_ == "push") {
        push(daughter, coord_func_->random_direction(wtl::prandom()));
    } else if (PACKING_ == "pushn") {
        push(daughter, to_nearest_empty(daughter->coord()));
    } else if (PACKING_ == "pushne") {
        pushn_everytime(daughter);
    } else if (PACKING_ == "fillpush") {
        fill_push(daughter, coord_func_->random_direction(wtl::prandom()));
    } else if (PACKING_ == "fill") {  //! @todo incorrect time scale
        return fill_empty(daughter);
    } else if (PACKING_ == "empty") {
        return insert_neighbor(daughter);
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
    std::shuffle(neighbors.begin(), neighbors.end(), wtl::prandom());
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
    daughter->set_coord(coord_func_->random_neighbor(daughter->coord(), wtl::prandom()));
    return tumor_.insert(daughter).second;
}

bool Tissue::swap_existing(std::shared_ptr<Cell>* x) {
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
    for (const auto& d: wtl::sample_knuth(coord_func_->directions(), search_max, wtl::prandom())) {
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

std::ostream& Tissue::write_segsites(std::ostream& ost, const std::vector<std::shared_ptr<Cell>>& subset) const {HERE;
    std::set<size_t> segsite_set;
    for (const auto p: subset) {
        const auto& sites = p->sites();
        segsite_set.insert(sites.begin(), sites.end());
    }
    const std::vector<size_t> segsites(segsite_set.begin(), segsite_set.end());
    std::vector<std::vector<size_t>> haplotypes;
    haplotypes.reserve(subset.size());
    for (const auto p: subset) {
        haplotypes.push_back(p->haplotype(segsites));
    }
    const size_t s = segsite_set.size();
    ost << "\n//\nsegsites: " << s << "\n";
    if (s > 0) {
        ost << "positions: ";
        wtl::ost_join(ost, std::vector<int>(s), " ") << "\n";
        for (const auto& x: haplotypes) {
            wtl::ost_join(ost, x, "") << "\n";
        }
    } else {ost << "\n";}
    return ost;
}

std::vector<std::shared_ptr<Cell>> Tissue::sample_random(const size_t n) const {HERE;
    std::vector<std::shared_ptr<Cell>> subset;
    if (tumor_.size() == stock_.size()) {  // death rate == 0
        subset = wtl::sample(stock_, n, wtl::prandom());
    } else {
        subset = wtl::sample(std::vector<std::shared_ptr<Cell>>(tumor_.begin(), tumor_.end()), n, wtl::prandom());
    }
    return subset;
}

std::vector<std::shared_ptr<Cell>>
Tissue::sample_if(std::function<bool(const std::vector<int>&)> predicate) const {HERE;
    std::vector<std::shared_ptr<Cell>> subset;
    for (const auto p: tumor_) {
        if (predicate(p->coord())) {subset.push_back(p);}
    }
    return subset;
}

std::string Tissue::snapshot_header() const {HERE;
    std::ostringstream oss;
    oss.precision(16);
    oss << "size" << sep_ << Cell::header(DIMENSIONS_, sep_);
    return oss.str();
}

std::string Tissue::snapshot() const {
    std::ostringstream oss;
    oss.precision(16);
    for (auto& item: stock_) {
        oss << tumor_.size() << sep_;
        item->write(oss, sep_);
    }
    return oss.str();
}

std::string Tissue::mutation_history() const {HERE;
    std::ostringstream oss;
    oss.precision(16);
    oss << "size" << sep_ << "mutant" << sep_ << "effect";
    std::vector<std::string> xyz{"x", "y", "z"};
    xyz.resize(DIMENSIONS_);
    for (auto x: xyz) {
        oss << sep_ << "origin_" << x;
    }
    oss << "\n";
    for (size_t i=0; i<mutation_coords_.size(); ++i) {
        oss << mutation_stages_[i] << sep_
            << Cell::MUTANT_IDS()[i] << sep_
            << Cell::MUTATION_EFFECTS()[i] << sep_;
        wtl::ost_join(oss, mutation_coords_[i], sep_) << "\n";
    }
    return oss.str();
}

std::string Tissue::evolution_history() const {
    return wtl::str_join(evolution_history_, "");
}

//! Stream operator for debug print
std::ostream& operator<< (std::ostream& ost, const Tissue& tissue) {
    return ost << tissue.tumor_;
}

template <class T> inline
void test_coordinate(const std::vector<int>& v) {HERE;
    T coord(v.size());
    std::cerr << coord.core() << std::endl;
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
    std::cerr << tissue.snapshot_header();
    std::cerr << tissue.snapshot() << std::endl;
    std::cerr << tissue.mutation_history() << std::endl;

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
