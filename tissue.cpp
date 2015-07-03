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

size_t Tissue::MAX_SIZE_ = 20000;
size_t Tissue::DIMENSIONS_ = 2;

boost::program_options::options_description& Tissue::opt_description() {
    namespace po = boost::program_options;
    static po::options_description desc{"Tissue"};
    desc.add_options()
        ("max,N", po::value<size_t>(&MAX_SIZE_)->default_value(MAX_SIZE_))
        ("dimensions,d", po::value<size_t>(&DIMENSIONS_)->default_value(DIMENSIONS_))
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
    tumor_.clear();
    while (tumor_.size() < n) {
        Gland daughter;
        daughter.mutate();
        std::vector<int> current_coords(DIMENSIONS_);
        push(std::move(daughter), &current_coords, random_direction(DIMENSIONS_));
        mutation_coords_.push_back(current_coords);
        mutation_stages_.push_back(tumor_.size());
    }
}

void Tissue::grow() {HERE;
    std::vector<std::vector<int> > coords;
    coords.reserve(MAX_SIZE_);
    for (const auto& item: tumor_) {
        coords.push_back(item.first);
    }
    while (tumor_.size() < MAX_SIZE_) {
        auto current_coords = *wtl::prandom().choice(coords.begin(), coords.end());
        auto& parent = tumor_[current_coords];
        Gland daughter = parent;
        if (parent.bernoulli_apoptosis()) {
            parent = std::move(daughter);
        } else {
            push(std::move(daughter), &current_coords, random_direction(DIMENSIONS_));
            coords.push_back(current_coords);
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
    axes.resize(DIMENSIONS_);
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
    xyz.resize(DIMENSIONS_);
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

void Tissue::unit_test() {
    std::cerr << __PRETTY_FUNCTION__ << std::endl;
    std::cerr.precision(15);
}
