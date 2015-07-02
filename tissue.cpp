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

size_t Tissue::MAX_SIZE_ = 10000;
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

std::string Tissue::grow() {HERE;
    std::ostringstream history;
    history.precision(16);
    history << "size\teffect\n";
    while (tumor_.size() < MAX_SIZE_) {
        auto parent = tumor_.begin();
        std::advance(parent, wtl::prandom().randrange(tumor_.size()));
        Gland daughter = parent->second;
        double fitness = daughter.fitness();
        if (daughter.mutate()) {
            history << tumor_.size() << "\t"
                    << daughter.fitness() - fitness << "\n";
        }
        if (parent->second.apoptosis()) {
            parent->second = daughter;
        } else {
            push(daughter, parent->first, random_direction(DIMENSIONS_));
        }
    }
    //derr(tumor_ << std::endl);
    return history.str();
}

void Tissue::push(const Gland& daughter, std::vector<int> coord, const std::vector<int>& direction) {
    std::transform(coord.begin(), coord.end(),direction.begin(), coord.begin(), std::plus<int>());
    auto result = tumor_.emplace(coord, daughter);
    if (!result.second) {
        push(result.first->second, coord, direction);
        tumor_[coord] = daughter;
    }
}

std::string Tissue::tsv() const {
    std::ostringstream ost;
    std::string sep("\t");
    ost.precision(16);
    std::vector<std::string> axes{"x", "y", "z"};
    for (size_t i=0; i<DIMENSIONS_; ++i) {
        ost << axes[i] << sep;
    }
    ost << "sites" << sep << "fitness\n";
    for (auto& item: tumor_) {
        wtl::ost_join(ost, item.first, sep) << sep;
        wtl::ost_join(ost, item.second.sites(), "|") << sep
            << item.second.fitness() << "\n";
    }
    return ost.str();
}

void Tissue::unit_test() {
    std::cerr << __PRETTY_FUNCTION__ << std::endl;
    std::cerr.precision(15);
}
