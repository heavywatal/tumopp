// -*- mode: c++; coding: utf-8 -*-
/*! @file tissue.cpp
    @brief Implementation of Gland class
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

boost::program_options::options_description& Tissue::opt_description() {
    namespace po = boost::program_options;
    static po::options_description desc{"Tissue"};
    desc.add_options()
    ;
    return desc;
}


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
    for (size_t i=1; i<n; ++i) {
        std::bitset<3> bs(i);
        std::vector<int> coord(dimensions);
        for (size_t j=0; j<n; ++j) {
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
            push(std::move(daughter), &current_coord, random_direction(current_coord));
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
    for (auto it=tumor_.begin(); tumor_.size() < max_size; ++it) {
        while (it != tumor_.end() && it->second.age() == age) {++it;}
        if (it == tumor_.end()) {
            if (tumor_.size() > max_size / 2) {
                return grow_random(max_size);
            }
            it = tumor_.begin();
            ++age;
        }
        it->second.stamp(age);
        auto current_coord = it->first;
        Gland daughter = it->second;
        if (it->second.bernoulli_apoptosis()) {
            it->second = std::move(daughter);
        } else {
            push(std::move(daughter), &current_coord, random_direction(current_coord));
        }
        if (Gland::bernoulli_mutation()) {
            tumor_[current_coord].mutate();
            mutation_coords_.push_back(std::move(current_coord));
            mutation_stages_.push_back(tumor_.size());
        }
    }
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

std::string Tissue::snapshot(const std::string& sep) const {HERE;
    std::ostringstream ost;
    ost.precision(16);
    std::vector<std::string> axes{"x", "y", "z"};
    axes.resize(coords_.front().size());
    wtl::ost_join(ost, axes, sep) << sep << "sites" << sep << "fitness\n";
    for (auto& item: tumor_) {
        wtl::ost_join(ost, item.first, sep) << sep;
        wtl::ost_join(ost, item.second.sites(), "|") << sep
            << item.second.fitness() << "\n";
    }
    return ost.str();
}

std::string Tissue::mutation_history(const std::string& sep) const {HERE;
    std::ostringstream ost;
    ost.precision(16);
    std::vector<std::string> xyz{"x", "y", "z"};
    xyz.resize(mutation_coords_.front().size());
    ost << "size" << sep << "effect";
    for (auto x: xyz) {
        ost << sep << "origin_" << x;
    }
    ost << "\n";
    for (size_t i=0; i<mutation_coords_.size(); ++i) {
        ost << mutation_stages_[i] << sep
            << Gland::MUTATION_EFFECTS()[i] << sep;
        wtl::ost_join(ost, mutation_coords_[i], sep) << "\n";
    }
    return ost.str();
}

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
}
