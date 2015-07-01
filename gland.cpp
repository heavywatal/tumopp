// -*- mode: c++; coding: utf-8 -*-
/*! @file gland.cpp
    @brief Implementation of Gland class
*/
#include "gland.h"

#include <cmath>
#include <iostream>
#include <sstream>

#include <boost/program_options.hpp>

#include "cxxwtils/prandom.hpp"
#include "cxxwtils/iostr.hpp"

double Gland::MUTATION_RATE_ = 1e-2;
double Gland::SIGMA_ = 1.0;
size_t Gland::MUTATED_SITES_ = 0;

boost::program_options::options_description& Gland::opt_description() {
    namespace po = boost::program_options;
    static po::options_description desc{"Gland"};
    desc.add_options()
        ("mutation,u", po::value<double>(&MUTATION_RATE_)->default_value(MUTATION_RATE_))
        ("sigma,s", po::value<double>(&SIGMA_)->default_value(SIGMA_))
    ;
    return desc;
}

void Gland::mutate() {
    if (wtl::prandom().bernoulli(MUTATION_RATE_)) {
        fitness_ += wtl::prandom().gauss(0.0, SIGMA_);
        fitness_ = std::max(fitness_, 1.0);
        fitness_ = std::min(fitness_, 5.0);
        sites_.push_back(MUTATED_SITES_++);
    }
}

bool Gland::apoptosis() const {
    return wtl::prandom().bernoulli(0.2 / fitness_);
}

int Gland::fission() const {
    if (malignant_ && !apoptosis()) {
        if (wtl::prandom().bernoulli(0.5)) {
            return 1;
        } else {
            return -1;
        }
    }
    return 0;
}


void Gland::unit_test() {
    std::cerr << __PRETTY_FUNCTION__ << std::endl;
    std::cerr.precision(15);
}
