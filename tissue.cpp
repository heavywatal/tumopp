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


size_t Tissue::DIMENSIONS_ = 2;
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


void Tissue::stain() {HERE;
    assert(tumor_.empty());
    for (const auto& coord: coord_func_->core()) {
        tumor_.insert(std::make_shared<Gland>(coord));
    }
}

void Tissue::grow(const size_t max_size) {HERE;
    assert(!tumor_.empty());
    evolution_history_.reserve(max_size);
    evolution_history_.push_back(snapshot());
    while (tumor_.size() < max_size) {
        std::vector<std::shared_ptr<Gland>> mothers;
        if (SCHEDULE_ == "random") {
            auto it = std::next(tumor_.begin(), wtl::prandom().randrange(tumor_.size()));
            mothers.push_back(*it);
        } else if (SCHEDULE_ == "poisson") {
            mothers.reserve(tumor_.size());
            for (auto it=tumor_.begin(); it!=tumor_.end(); ++it) {
                //! @todo fitness
                const size_t n = wtl::prandom().poisson(1.0);
                for (size_t i=0; i<n; ++i) {
                    mothers.push_back(*it);
                }
            }
        } else if (SCHEDULE_ == "even") {
            mothers.reserve(tumor_.size());
            for (auto it=tumor_.begin(); it!=tumor_.end(); ++it) {
                mothers.push_back(*it);
            }
        } else {exit(1);}
        for (auto& mother: mothers) {
            if (mother->bernoulli_birth()) {
                auto daughter = std::make_shared<Gland>(*mother);  // Gland copy ctor
                bool success = true;
                if (PACKING_ == "push") {
                    push(daughter, coord_func_->random_direction(wtl::prandom()));
                } else if (PACKING_ == "fillpush") {
                    fill_push(daughter, coord_func_->random_direction(wtl::prandom()));
                } else if (PACKING_ == "fillwalk") {
                    fill_walk(daughter);
                } else if (PACKING_ == "fill") {
                    success = fill_empty(daughter);
                } else if (PACKING_ == "empty") {
                    success = insert_neighbor(daughter);
                }
                if (!success) {continue;}
                if (Gland::bernoulli_mutation()) {
                    daughter->mutate();
                    mutation_coords_.push_back(daughter->coord());
                    mutation_stages_.push_back(tumor_.size());
                }
            }
            if (mother->bernoulli_death()) {
                tumor_.erase(mother);
            } else if (Gland::bernoulli_mutation()) {
                mother->mutate();
                mutation_coords_.push_back(mother->coord());
                mutation_stages_.push_back(tumor_.size());
            }
            if (tumor_.size() <= 128) {
                evolution_history_.push_back(snapshot());
            }
        }
    }
}

std::string Tissue::evolution_history() const {
    return wtl::str_join(evolution_history_, "");
}

void Tissue::push(const std::shared_ptr<Gland> moving, const std::vector<int>& direction) {
    moving->set_coord(moving->coord() + direction);
    const auto result = tumor_.insert(moving);
    if (!result.second) {
        const std::shared_ptr<Gland> existing = *result.first;
        tumor_.erase(result.first);
        assert(tumor_.insert(moving).second);
        push(existing, direction);
    }
}

std::shared_ptr<Gland> Tissue::push_pop(const std::shared_ptr<Gland> x) {
    const auto pos = tumor_.find(x);
    assert(pos != tumor_.end());
    const std::shared_ptr<Gland> existing = *pos;
    tumor_.erase(pos);
    assert(tumor_.insert(x).second);
    return existing;
}

bool Tissue::fill_empty(const std::shared_ptr<Gland> moving) {
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

void Tissue::fill_push(std::shared_ptr<Gland> moving, const std::vector<int>& direction) {
    while (!fill_empty(moving)) {
        moving->set_coord(moving->coord() + direction);
        moving = push_pop(moving);
    }
}

void Tissue::fill_walk(std::shared_ptr<Gland> moving) {
    while (!fill_empty(moving)) {
        moving->set_coord(coord_func_->random_neighbor(moving->coord(), wtl::prandom()));
        moving = push_pop(moving);
    }
}

bool Tissue::insert_neighbor(const std::shared_ptr<Gland> daughter) {
    daughter->set_coord(coord_func_->random_neighbor(daughter->coord(), wtl::prandom()));
    const auto result = tumor_.insert(daughter);
    return result.second;
}

std::string Tissue::snapshot_header() const {HERE;
    std::ostringstream ost;
    ost.precision(16);
    ost << "time" << sep_ << "size" << sep_ << Gland::header(DIMENSIONS_, sep_);
    return ost.str();
}

std::string Tissue::snapshot() const {
    std::ostringstream ost;
    ost.precision(16);
    for (auto& item: tumor_) {
        ost << evolution_history_.size() << sep_ << tumor_.size() << sep_;
        item->write(ost, sep_);
    }
    return ost.str();
}

std::string Tissue::mutation_history() const {HERE;
    std::ostringstream ost;
    ost.precision(16);
    ost << "size" << sep_ << "mutant" << sep_ << "effect";
    std::vector<std::string> xyz{"x", "y", "z"};
    xyz.resize(DIMENSIONS_);
    for (auto x: xyz) {
        ost << sep_ << "origin_" << x;
    }
    ost << "\n";
    for (size_t i=0; i<mutation_coords_.size(); ++i) {
        ost << mutation_stages_[i] << sep_
            << Gland::MUTANT_IDS()[i] << sep_
            << Gland::MUTATION_EFFECTS()[i] << sep_;
        wtl::ost_join(ost, mutation_coords_[i], sep_) << "\n";
    }
    return ost.str();
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
    tissue.stain();
    tissue.grow(100);
    std::cerr << tissue.size() << ": " << tissue.radius() << std::endl;
    tissue.grow(1000);
    std::cerr << tissue.size() << ": " << tissue.radius() << std::endl;
}

void Tissue::unit_test() {HERE;
    std::cerr.precision(15);

    Tissue tissue;
    tissue.stain();
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
