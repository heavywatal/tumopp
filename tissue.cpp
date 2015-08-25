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

//! Program options
boost::program_options::options_description& Tissue::opt_description() {
    namespace po = boost::program_options;
    static po::options_description desc{"Tissue"};
    desc.add_options()
    ;
    return desc;
}

template <class T>
std::vector<T> operator+(const std::vector<T>& lhs, const std::vector<T>& rhs) {
    assert(lhs.size() == rhs.size());
    std::vector<T> result;
    result.reserve(lhs.size());
    std::transform(lhs.begin(), lhs.end(), rhs.begin(),
                   std::back_inserter(result), std::plus<T>());
    return result;
}

template <class T>
std::vector<T>& operator+=(std::vector<T>& lhs, const std::vector<T>& rhs) {
    assert(lhs.size() == rhs.size());
    std::transform(lhs.begin(), lhs.end(), rhs.begin(),
                   lhs.begin(), std::plus<T>());
    return lhs;
}

//! @arg dimensions 2 or 3
inline std::vector<std::vector<int>> all_directions(const size_t dimensions) {
    assert(dimensions == 2 || dimensions == 3);
    std::vector<std::vector<int>> output;
    output.reserve(std::pow(3, dimensions) - 1);
    for (const int x: {-1, 0, 1}) {
        for (const int y: {-1, 0, 1}) {
            if (dimensions == 2) {
                if (x == 0 && y == 0) continue;
                output.push_back({x, y});
                continue;
            }
            for (const int z: {-1, 0, 1}) {
                if (x == 0 && y == 0 && z == 0) continue;
                output.push_back({x, y, z});
            }
        }
    }
    return output;
}

//! @arg dimensions 2 or 3
inline std::vector<std::vector<int>> proximal_directions(const size_t dimensions) {
    assert(dimensions == 2 || dimensions == 3);
    std::vector<std::vector<int>> output;
    output.reserve(2 * dimensions);
    std::vector<int> v(dimensions, 0);
    v.back() += 1;
    do {
        output.push_back(v);
    } while (std::next_permutation(v.begin(), v.end()));
    v.assign(dimensions, 0);
    v.front() -= 1;
    do {
        output.push_back(v);
    } while (std::next_permutation(v.begin(), v.end()));
    return output;
}

//! @arg dimensions 2
inline std::vector<std::vector<int>> hex_directions(const size_t dimensions) {
    assert(dimensions == 2);
    std::vector<std::vector<int>> output;
    output.reserve(std::pow(3, dimensions) - 1);
    for (const int x: {-1, 0, 1}) {
        for (const int y: {-1, 0, 1}) {
            if (dimensions == 2) {
                if (x == 0 && y == 0) continue;
                output.push_back({x, y});
                continue;
            }
            for (const int z: {-1, 0, 1}) {
                if (x == 0 && y == 0 && z == 0) continue;
                output.push_back({x, y, z});
            }
        }
    }
    return output;
}

//! random vector of {-1, 0, 1}
inline std::vector<int> random_direction(const std::vector<int>& current_coord) {
    const size_t dimensions = current_coord.size();
    std::vector<int> direction(dimensions, -1);
    for (size_t i=0; i<dimensions; ++i) {
        direction[i] += wtl::prandom().randrange(3);
    }
    if (direction == std::vector<int>(dimensions, 0)) {
        return random_direction(current_coord);
    }
    return direction;
}

//! random vector of {-1, 0, 1} outward
inline std::vector<int> random_outward(const std::vector<int>& current_coord) {
    const size_t dimensions = current_coord.size();
    std::vector<int> direction(dimensions);
    for (size_t i=0; i<dimensions; ++i) {
        int n = wtl::prandom().randrange(2);
        if (current_coord[i] > 0) {direction[i] = n;}
        else {direction[i] = -n;}
    }
    if (direction == std::vector<int>(dimensions)) {
        return random_outward(current_coord);
    }
    return direction;
}

void Tissue::init_regularly() {HERE;
    const size_t dimensions = coords_.front().size();
    const size_t n = std::pow(2, dimensions);
    coords_.reserve(n);
    for (size_t i=tumor_.size(); i<n; ++i) {
        std::bitset<3> bs(i);
        std::vector<int> coord(dimensions);
        for (size_t j=0; j<dimensions; ++j) {
            coord[j] = static_cast<int>(bs[j]);
        }
        emplace(coord, Gland());
    }
}

void Tissue::init_randomly() {HERE;
    const size_t dimensions = coords_.front().size();
    const size_t n = std::pow(2, dimensions);
    while (tumor_.size() < n) {
        std::vector<int> coord = coords_.front();
        push(Gland(), &coord, random_direction(coord));
    }
}

void Tissue::stain() {HERE;
    for (auto& item: tumor_) {
        item.second.mutate();
        mutation_coords_.push_back(item.first);
        mutation_stages_.push_back(mutation_coords_.size());
    }
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
            fill_neighbor(std::move(daughter), current_coord);
//            push(std::move(daughter), &current_coord, random_direction(current_coord));
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
            push_neighbor(std::move(daughter), current_coord);
//            push(std::move(daughter), &current_coord, random_direction(current_coord));
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
    std::transform(coord->begin(), coord->end(), direction.begin(), coord->begin(), std::plus<int>());
    auto it = tumor_.find(*coord);
    if (it == tumor_.end()) {
        emplace(*coord, std::move(daughter));
    } else {
        push(std::move(it->second), coord, direction);
        it->second = std::move(daughter);
        *coord = it->first;
    }
}

//! @todo
void Tissue::hex_neighbor(Gland&& daughter, const std::vector<int>& current_coord,
                           const std::vector<int>& direction) {
    const size_t dimensions = current_coord.size();
//    static const auto directions = all_directions(dimensions);
//    static const auto directions = proximal_directions(dimensions);
    static const auto directions = Hex::directions();
//    auto neighbors = directions;
    auto neighbors = Hex(current_coord).neighbors();
    std::vector<std::vector<int>> empty_neighbors;
    empty_neighbors.reserve(neighbors.size());
    for (auto& x: neighbors) {
//        x += current_coord;
        if (tumor_.find(x.vec()) == tumor_.end()) {
            empty_neighbors.push_back(x.vec());
        }
    }
    std::vector<int> new_coord;
    if (empty_neighbors.empty()) {
        auto new_dir = direction;
        if (direction.empty()) {
            new_dir = wtl::prandom().choice(directions.begin(), directions.end())->vec();
        }
        new_coord = current_coord + new_dir;
        push_neighbor(std::move(tumor_[current_coord]), new_coord, new_dir);
    } else {
        new_coord = *std::min_element(empty_neighbors.begin(), empty_neighbors.end(),
            [](const std::vector<int>& x, const std::vector<int>& y){
                return wtl::devsq(x) < wtl::devsq(y);
        });
        emplace(new_coord, std::move(tumor_[current_coord]));
    }
    tumor_[current_coord] = std::move(daughter);
}

//! @todo
void Tissue::push_neighbor(Gland&& daughter, const std::vector<int>& current_coord,
                           const std::vector<int>& direction) {
    const size_t dimensions = current_coord.size();
    static const auto directions = all_directions(dimensions);
//    static const auto directions = proximal_directions(dimensions);
    auto neighbors = directions;
    std::vector<std::vector<int>> empty_neighbors;
    empty_neighbors.reserve(neighbors.size());
    for (auto& x: neighbors) {
        x += current_coord;
        if (tumor_.find(x) == tumor_.end()) {
            empty_neighbors.push_back(x);
        }
    }
    std::vector<int> new_coord;
    if (empty_neighbors.empty()) {
        auto new_dir = direction;
        if (direction.empty()) {
            new_dir = *wtl::prandom().choice(directions.begin(), directions.end());
        }
        new_coord = current_coord + new_dir;
        push_neighbor(std::move(tumor_[current_coord]), new_coord, new_dir);
    } else {
        new_coord = *std::min_element(empty_neighbors.begin(), empty_neighbors.end(),
            [](const std::vector<int>& x, const std::vector<int>& y){
                return wtl::devsq(x) < wtl::devsq(y);
        });
        emplace(new_coord, std::move(tumor_[current_coord]));
    }
    tumor_[current_coord] = std::move(daughter);
}

//! @todo
void Tissue::fill_neighbor(Gland&& daughter, const std::vector<int>& current_coord) {
    const size_t dimensions = current_coord.size();
    static const auto directions = all_directions(dimensions);
//    static const auto directions = proximal_directions(dimensions);
    auto neighbors = directions;
    std::vector<std::vector<int>> empty_neighbors;
    empty_neighbors.reserve(neighbors.size());
    for (auto& x: neighbors) {
        x += current_coord;
        if (tumor_.find(x) == tumor_.end()) {
            empty_neighbors.push_back(x);
        }
    }
    std::vector<int> new_coord;
    if (empty_neighbors.empty()) {
        new_coord = *wtl::prandom().choice(neighbors.begin(), neighbors.end());
        fill_neighbor(std::move(tumor_[current_coord]), new_coord);
    } else {
        new_coord = *std::min_element(empty_neighbors.begin(), empty_neighbors.end(),
            [](const std::vector<int>& x, const std::vector<int>& y){
                return wtl::devsq(x) < wtl::devsq(y);
        });
        emplace(new_coord, std::move(tumor_[current_coord]));
    }
    tumor_[current_coord] = std::move(daughter);
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

void Tissue::unit_test() {
    std::cerr << __PRETTY_FUNCTION__ << std::endl;
    std::cerr.precision(15);
    Tissue tissue(3);
    tissue.stain();
    tissue.grow_even(10);
    std::cerr << tissue << std::endl;
    std::cerr << tissue.coords_ << std::endl;
    std::cerr << tissue.snapshot() << std::endl;
    std::cerr << tissue.mutation_history() << std::endl;
    std::cerr << proximal_directions(2) << std::endl;
    std::cerr << proximal_directions(3) << std::endl;
    Hex h(3, 2);
    std::cerr << h << std::endl;
    std::cerr << h.radius() << std::endl;
    std::cerr << Hex::distance(h, Hex(1, 3)) << std::endl;
    std::cerr << Hex::directions() << std::endl;
    std::cerr << h.neighbors() << std::endl;
}
