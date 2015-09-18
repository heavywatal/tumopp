// -*- mode: c++; coding: utf-8 -*-
/*! @file tissue.cpp
    @brief Implementation of Tissue class
*/
#include "tissue.h"

#include <cassert>
#include <iostream>
#include <sstream>

#include <boost/program_options.hpp>

#include "cxxwtils/prandom.hpp"
#include "cxxwtils/iostr.hpp"
#include "cxxwtils/gz.hpp"
#include "cxxwtils/debug.hpp"
#include "cxxwtils/algorithm.hpp"


size_t Tissue::DIMENSIONS_ = 2;

//! Program options
/*! @return Program options description

    Command line option | Symbol | Variable
    --------------------| ------ | -------------------------
    `-D,--dimensions`   | -      | Tissue::DIMENSIONS_
*/
boost::program_options::options_description& Tissue::opt_description() {
    namespace po = boost::program_options;
    static po::options_description desc{"Tissue"};
    desc.add_options()
        ("dimensions,D", po::value<size_t>(&DIMENSIONS_)->default_value(DIMENSIONS_))
    ;
    return desc;
}

void Tissue::grow_random(const size_t max_size) {HERE;
    evolution_history_.reserve(max_size);
    evolution_history_.push_back(snapshot());
    while (tumor_.size() < max_size) {
        auto it = std::next(tumor_.begin(), wtl::prandom().randrange(tumor_.size()));
        auto& mother = *it;
        auto daughter = std::make_shared<Gland>(*mother);
        if (mother->bernoulli_apoptosis()) {
            tumor_.erase(it);
            tumor_.insert(daughter);
        } else {
            //push(daughter);
            push_fill(daughter);
            //walk_fill(daughter);
        }
        if (Gland::bernoulli_mutation()) {
            daughter->mutate();
            mutation_coords_.push_back(daughter->coord());
            mutation_stages_.push_back(tumor_.size());
        }
        if (tumor_.size() <= 128) {
            evolution_history_.push_back(snapshot());
        }
    }
}

void Tissue::grow_even(const size_t max_size) {HERE;
    evolution_history_.reserve(max_size);
    evolution_history_.push_back(snapshot());
    while (tumor_.size() < max_size) {
        std::vector<std::shared_ptr<Gland>> new_glands;
        new_glands.reserve(tumor_.size());
        for (auto it=tumor_.begin(); it!=tumor_.end(); ++it) {
            auto& mother = *it;
            if (mother->bernoulli_apoptosis()) {
                if (Gland::bernoulli_mutation()) {
                    mother->mutate();
                    mutation_coords_.push_back(mother->coord());
                    mutation_stages_.push_back(tumor_.size());
                }
            } else {
                new_glands.push_back(mother);
            }
        }
        for (auto& gland: new_glands) {
            auto daughter = std::make_shared<Gland>(*gland);  // Gland copy ctor
            if (Gland::bernoulli_mutation()) {
                daughter->mutate();
                mutation_coords_.push_back(daughter->coord());
                mutation_stages_.push_back(tumor_.size());
            }
            //push(daughter);
            push_fill(daughter);
            //walk_fill(daughter);
            if (tumor_.size() <= 128) {
                evolution_history_.push_back(snapshot());
            }
        }
    }
}

void Tissue::grow_poisson(const size_t max_size) {HERE;
    //! @todo
}

std::string Tissue::evolution_history() const {
    return wtl::str_join(evolution_history_, "");
}


void Tissue::push(const std::shared_ptr<Gland>& daughter, const std::vector<int>& direction) {
    auto new_dir = direction;
    if (direction.empty()) {
        const auto directions = coord_func_->directions(daughter->coord().size());
        new_dir = *wtl::prandom().choice(directions.begin(), directions.end());
    }
    auto new_coord = daughter->coord() + new_dir;
    auto found = tumor_.find(daughter);
    auto existing = *found;
    tumor_.erase(found);
    existing->set_coord(new_coord);
    if (!tumor_.insert(existing).second) {
        push(existing, new_dir);
    }
    tumor_.insert(daughter);
}


void Tissue::push_fill(const std::shared_ptr<Gland>& daughter, const std::vector<int>& direction) {
    static const auto directions = coord_func_->directions(daughter->coord().size());
    auto neighbors = coord_func_->neighbors(daughter->coord());
    if (!direction.empty()) {
        auto previous = std::find(neighbors.begin(), neighbors.end(), daughter->coord() - direction);
        neighbors.erase(previous);
    }
    std::vector<std::vector<int>> empty_neighbors;
    empty_neighbors.reserve(neighbors.size());
    for (auto& x: neighbors) {
        if (tumor_.find(std::shared_ptr<Gland>(new Gland(x))) == tumor_.end()) {
            empty_neighbors.push_back(x);
        }
    }
    auto found = tumor_.find(daughter);
    auto existing = *found;
    tumor_.erase(found);
    if (empty_neighbors.empty()) {
        auto new_dir = direction;
        if (direction.empty()) {
            new_dir = *wtl::prandom().choice(directions.begin(), directions.end());
        }
        auto new_coord = daughter->coord() + new_dir;
        existing->set_coord(new_coord);
        push_fill(existing, new_dir);
    } else {
        auto new_coord = *std::min_element(empty_neighbors.begin(), empty_neighbors.end(),
            [this](const std::vector<int>& x, const std::vector<int>& y) {
                return coord_func_->distance(x) < coord_func_->distance(y);
        });
        existing->set_coord(new_coord);
        tumor_.insert(existing);
    }
    tumor_.insert(daughter);
}

void Tissue::walk_fill(const std::shared_ptr<Gland>& daughter) {
    static const auto directions = coord_func_->directions(daughter->coord().size());
    auto neighbors = coord_func_->neighbors(daughter->coord());
    std::vector<std::vector<int>> empty_neighbors;
    empty_neighbors.reserve(neighbors.size());
    for (auto& x: neighbors) {
        if (tumor_.find(std::shared_ptr<Gland>(new Gland(x))) == tumor_.end()) {
            empty_neighbors.push_back(x);
        }
    }
    auto found = tumor_.find(daughter);
    auto existing = *found;
    tumor_.erase(found);
    if (empty_neighbors.empty()) {
        auto new_coord = *wtl::prandom().choice(neighbors.begin(), neighbors.end());
        existing->set_coord(new_coord);
        walk_fill(existing);
    } else {
        auto new_coord = *std::min_element(empty_neighbors.begin(), empty_neighbors.end(),
            [this](const std::vector<int>& x, const std::vector<int>& y) {
                return coord_func_->distance(x) < coord_func_->distance(y);
        });
        existing->set_coord(new_coord);
        tumor_.insert(existing);
    }
    tumor_.insert(daughter);
}

//! @todo
void Tissue::push_layer(const std::shared_ptr<Gland>& daughter) {
}


void Tissue::stain() {HERE;
    assert(tumor_.empty());
    for (const auto& coord: coord_func_->origins(DIMENSIONS_)) {
        std::shared_ptr<Gland> founder(new Gland(coord));
        founder->mutate();
        tumor_.insert(founder);
        mutation_coords_.push_back(coord);
        mutation_stages_.push_back(mutation_coords_.size());
    }
}

std::string Tissue::snapshot_header() const {HERE;
    std::ostringstream ost;
    ost.precision(16);
    ost << "time" << sep_ << "size" << sep_;
    std::vector<std::string> axes{"x", "y", "z"};
    axes.resize(DIMENSIONS_);
    wtl::ost_join(ost, axes, sep_) << sep_ << "sites" << sep_ << "fitness\n";
    return ost.str();
}

std::string Tissue::snapshot() const {
    std::ostringstream ost;
    ost.precision(16);
    for (auto& item: tumor_) {
        ost << evolution_history_.size() << sep_ << tumor_.size() << sep_;
        wtl::ost_join(ost, item->coord(), sep_) << sep_;
        wtl::ost_join(ost, item->sites(), "|") << sep_
            << item->fitness() << "\n";
    }
    return ost.str();
}

std::string Tissue::mutation_history() const {HERE;
    std::ostringstream ost;
    ost.precision(16);
    std::vector<std::string> xyz{"x", "y", "z"};
    xyz.resize(mutation_coords_.front().size());
    ost << "size" << sep_ << "effect";
    for (auto x: xyz) {
        ost << sep_ << "origin_" << x;
    }
    ost << "\n";
    for (size_t i=0; i<mutation_coords_.size(); ++i) {
        ost << mutation_stages_[i] << sep_
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
void test_coordinate(const std::vector<int>& v) {
    T coord;
    std::cerr << typeid(coord).name() << std::endl;
    std::cerr << coord.distance(v) << std::endl;
    std::cerr << coord.neighbors(v) << std::endl;
    std::cerr << coord.directions(v.size()) << std::endl;
}

void Tissue::unit_test() {
    std::cerr << __PRETTY_FUNCTION__ << std::endl;
    std::cerr.precision(15);

    test_coordinate<Lattice>({3, -2});
    test_coordinate<Diagonal>({3, -2});
    test_coordinate<Hexagonal>({3, -2});

    test_coordinate<Lattice>({3, -2, 1});
    test_coordinate<Diagonal>({3, -2, 1});
    test_coordinate<Hexagonal>({3, -2, 1});

    std::cerr << Lattice().origins(2) << std::endl;
    std::cerr << Lattice().origins(3) << std::endl;
    std::cerr << Hexagonal().origins(2) << std::endl;
    std::cerr << Hexagonal().origins(3) << std::endl;

    Tissue tissue;
    tissue.stain();
    tissue.grow_even(10);
    std::cerr << tissue << std::endl;
    std::cerr << tissue.snapshot_header();
    std::cerr << tissue.snapshot() << std::endl;
    std::cerr << tissue.mutation_history() << std::endl;
}
