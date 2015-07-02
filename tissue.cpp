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

size_t Tissue::MAX_SIZE_ = 100;
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
    std::vector<int> direction(dimension);
    int sign = wtl::prandom().randrange(2);  // 0 or 1
    sign *= 2;  // 0 or 2
    sign -= 1;  // -1 or 1
    direction[wtl::prandom().randrange(dimension)] = sign;
    return direction;
}

void Tissue::grow() {
    tumor_.emplace(std::vector<int>(DIMENSIONS_), Gland());
    std::cout << tumor_ << std::endl;
    std::ostringstream history;
    history.precision(16);
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
    std::cout << tumor_ << std::endl;
    std::cout << history.str() << std::endl;
}

void Tissue::push(const Gland& daughter, std::vector<int> coord, const std::vector<int>& direction) {
    std::transform(coord.begin(), coord.end(),direction.begin(), coord.begin(), std::plus<int>());
    auto result = tumor_.emplace(coord, daughter);
    if (!result.second) {
        push(result.first->second, coord, direction);
        tumor_[coord] = daughter;
    }
}

void Tissue::unit_test() {
    std::cerr << __PRETTY_FUNCTION__ << std::endl;
    std::cerr.precision(15);
}
