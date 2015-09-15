// -*- mode: c++; coding: utf-8 -*-
/*! @file tissue.cpp
    @brief Implementation of Tissue class
*/
#include "tissue.h"

#include <cmath>
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

Tissue::Tissue(const std::vector<int>& origin): coord_func_(new Hex()) {
    emplace(origin, Gland());
    init_regularly();
}

void Tissue::grow_random(const size_t max_size) {HERE;
    coords_.reserve(max_size);
    while (tumor_.size() < max_size) {
        auto current_coord = *wtl::prandom().choice(coords_.begin(), coords_.end());
        auto& parent = tumor_[current_coord];
        Gland daughter = parent;
        if (parent.bernoulli_apoptosis()) {
            parent = std::move(daughter);
        } else {
//            push(std::move(daughter), &current_coord);
            push_fill(std::move(daughter), current_coord);
//            walk_fill(std::move(daughter), current_coord);
        }
        if (Gland::bernoulli_mutation()) {
            tumor_[current_coord].mutate();
            mutation_coords_.push_back(std::move(current_coord));
            mutation_stages_.push_back(tumor_.size());
        }
    }
}

void Tissue::grow_even(const size_t max_size) {HERE;
    coords_.reserve(max_size);
    size_t age = 1;
    evolution_history_.reserve(max_size);
    evolution_history_.push_back(snapshot());
    for (auto it=coords_.rbegin(); tumor_.size() < max_size; ++it) {
        while (it != coords_.rend() && tumor_[*it].age() == age) {++it;}
        if (it == coords_.rend()) {
            it = coords_.rbegin();
            ++age;
        }
        auto& mother = tumor_[*it];
        mother.stamp(age);
        auto current_coord = *it;
        Gland daughter = mother;
        if (mother.bernoulli_apoptosis()) {
            mother = std::move(daughter);
        } else {
//            push(std::move(daughter), &current_coord);
            push_fill(std::move(daughter), current_coord);
//            walk_fill(std::move(daughter), current_coord);
        }
        if (Gland::bernoulli_mutation()) {
            tumor_[current_coord].mutate();
            mutation_coords_.push_back(std::move(current_coord));
            mutation_stages_.push_back(tumor_.size());
        }
        if (tumor_.size() <= 256) {
            evolution_history_.push_back(snapshot());
        }
    }
}

void Tissue::grow_poisson(const size_t max_size) {HERE;
    //! @todo
}

std::string Tissue::evolution_history() const {
    return wtl::str_join(evolution_history_, "");
}

void Tissue::emplace(const std::vector<int>& coord, Gland&& daughter) {
    coords_.push_back(coord);
    tumor_.emplace(coord, std::move(daughter));
}

void Tissue::push(Gland&& daughter, std::vector<int>* coord, const std::vector<int>& direction) {
    // TODO: rate-limiting
    auto new_dir = direction;
    if (direction.empty()) {
        const auto directions = coord_func_->directions(coord->size());
        new_dir = *wtl::prandom().choice(directions.begin(), directions.end());
    }
    *coord += new_dir;
    auto it = tumor_.find(*coord);
    if (it == tumor_.end()) {
        emplace(*coord, std::move(daughter));
    } else {
        push(std::move(it->second), coord, new_dir);
        it->second = std::move(daughter);
        *coord = it->first;
    }
}

//! @todo
void Tissue::push_fill(Gland&& daughter, const std::vector<int>& current_coord,
                       const std::vector<int>& direction) {
    static const auto directions = coord_func_->directions(current_coord.size());
    const auto neighbors = coord_func_->neighbors(current_coord);
    std::vector<std::vector<int>> empty_neighbors;
    empty_neighbors.reserve(neighbors.size());
    for (auto& x: neighbors) {
        if (tumor_.find(x) == tumor_.end()) {
            empty_neighbors.push_back(x);
        }
    }
    if (empty_neighbors.empty()) {
        auto new_dir = direction;
        if (direction.empty()) {
            new_dir = *wtl::prandom().choice(directions.begin(), directions.end());
        }
        auto new_coord = current_coord + new_dir;
        push_fill(std::move(tumor_[current_coord]), new_coord, new_dir);
    } else {
        auto new_coord = *std::min_element(empty_neighbors.begin(), empty_neighbors.end(),
            [this](const std::vector<int>& x, const std::vector<int>& y) {
                return coord_func_->distance(x) < coord_func_->distance(y);
        });
        emplace(new_coord, std::move(tumor_[current_coord]));
    }
    tumor_[current_coord] = std::move(daughter);
}

//! @todo
void Tissue::walk_fill(Gland&& daughter, const std::vector<int>& current_coord) {
    static const auto directions = coord_func_->directions(current_coord.size());
    const auto neighbors = coord_func_->neighbors(current_coord);
    std::vector<std::vector<int>> empty_neighbors;
    empty_neighbors.reserve(neighbors.size());
    for (auto& x: neighbors) {
        if (tumor_.find(x) == tumor_.end()) {
            empty_neighbors.push_back(x);
        }
    }
    if (empty_neighbors.empty()) {
        auto new_coord = *wtl::prandom().choice(neighbors.begin(), neighbors.end());
        walk_fill(std::move(tumor_[current_coord]), new_coord);
    } else {
        auto new_coord = *std::min_element(empty_neighbors.begin(), empty_neighbors.end(),
            [this](const std::vector<int>& x, const std::vector<int>& y) {
                return coord_func_->distance(x) < coord_func_->distance(y);
        });
        emplace(new_coord, std::move(tumor_[current_coord]));
    }
    tumor_[current_coord] = std::move(daughter);
}

void Tissue::init_regularly() {HERE;
    for (const auto& coord: coord_func_->neighbors(coords_.front())) {
        emplace(coord, Gland());
    }
}

void Tissue::stain() {HERE;
    for (auto& item: tumor_) {
        item.second.mutate();
        mutation_coords_.push_back(item.first);
        mutation_stages_.push_back(mutation_coords_.size());
    }
}

std::string Tissue::snapshot_header() const {HERE;
    std::ostringstream ost;
    ost.precision(16);
    ost << "time" << sep_ << "size" << sep_;
    std::vector<std::string> axes{"x", "y", "z"};
    axes.resize(coords_.front().size());
    wtl::ost_join(ost, axes, sep_) << sep_ << "sites" << sep_ << "fitness\n";
    return ost.str();
}

std::string Tissue::snapshot() const {
    std::ostringstream ost;
    ost.precision(16);
    for (auto& item: tumor_) {
        ost << evolution_history_.size() << sep_ << tumor_.size() << sep_;
        wtl::ost_join(ost, item.first, sep_) << sep_;
        wtl::ost_join(ost, item.second.sites(), "|") << sep_
            << item.second.fitness() << "\n";
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
    Tissue tissue;
    tissue.stain();
    tissue.grow_even(10);
    std::cerr << tissue << std::endl;
    std::cerr << tissue.coords_ << std::endl;
    std::cerr << tissue.snapshot() << std::endl;
    std::cerr << tissue.mutation_history() << std::endl;

    test_coordinate<Orthogonal>({3, -2});
    test_coordinate<Lattice>({3, -2});
    test_coordinate<Hex>({3, -2});

    test_coordinate<Orthogonal>({3, -2, 1});
    test_coordinate<Lattice>({3, -2, 1});
    test_coordinate<Hex>({3, -2, 1});
}
