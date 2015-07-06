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


inline std::vector<int> random_direction(const size_t dimension) {
    std::vector<int> direction(dimension, -1);
    for (size_t i=0; i<dimension; ++i) {
        direction[i] += wtl::prandom().randrange(3);
    }
    if (direction == std::vector<int>(dimension, 0)) {
        return random_direction(dimension);
    }
    return direction;
}

void Tissue::mark(const size_t n) {HERE;
    while (tumor_.size() < n) {
        std::vector<int> coord = coords_.front();
        push(Gland(), &coord, random_direction(coord.size()));
    }
    for (auto& item: tumor_) {
        item.second.mutate();
        mutation_coords_.push_back(item.first);
        mutation_stages_.push_back(mutation_coords_.size());
    }
}

void Tissue::grow(const size_t max_size) {HERE;
    coords_.reserve(max_size);
    while (tumor_.size() < max_size) {
        auto current_coords = *wtl::prandom().choice(coords_.begin(), coords_.end());
        auto& parent = tumor_[current_coords];
        Gland daughter = parent;
        if (parent.bernoulli_apoptosis()) {
            parent = std::move(daughter);
        } else {
            push(std::move(daughter), &current_coords, random_direction(current_coords.size()));
        }
        if (Gland::bernoulli_mutation()) {
            tumor_[current_coords].mutate();
            mutation_coords_.push_back(std::move(current_coords));
            mutation_stages_.push_back(tumor_.size());
        }
    }
}

void Tissue::push(Gland&& daughter, std::vector<int>* coord, const std::vector<int>& direction) {
    // TODO: rate-limiting
    std::transform(coord->begin(), coord->end(), direction.begin(), coord->begin(), std::plus<int>());
    auto it = tumor_.find(*coord);
    if (it == tumor_.end()) {
        coords_.push_back(*coord);
        tumor_.emplace(*coord, std::move(daughter));
    } else {
        push(std::move(it->second), coord, direction);
        it->second = std::move(daughter);
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
    tissue.mark(4);
    tissue.grow(10);
    std::cerr << tissue << std::endl;
    std::cerr << tissue.coords_ << std::endl;
    std::cerr << tissue.snapshot() << std::endl;
    std::cerr << tissue.mutation_history() << std::endl;
}
