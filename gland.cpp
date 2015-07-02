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
double Gland::MUTATION_SIGMA_ = 1.0;
double Gland::APOPTOSIS_RATE_ = 0.2;

size_t Gland::MUTATED_SITES_ = 0;

boost::program_options::options_description& Gland::opt_description() {
    namespace po = boost::program_options;
    static po::options_description desc{"Gland"};
    desc.add_options()
        ("mutation,u", po::value<double>(&MUTATION_RATE_)->default_value(MUTATION_RATE_))
        ("sigma,s", po::value<double>(&MUTATION_SIGMA_)->default_value(MUTATION_SIGMA_))
        ("apoptosis,a", po::value<double>(&APOPTOSIS_RATE_)->default_value(APOPTOSIS_RATE_))
    ;
    return desc;
}

bool Gland::mutate() {
    if (wtl::prandom().bernoulli(MUTATION_RATE_)) {
        fitness_ += wtl::prandom().gauss(0.0, MUTATION_SIGMA_);
        fitness_ = std::max(fitness_, 1.0);
        fitness_ = std::min(fitness_, 5.0);
        sites_.push_back(MUTATED_SITES_++);
        return true;
    }
    return false;
}

bool Gland::apoptosis() const {
    return wtl::prandom().bernoulli(APOPTOSIS_RATE_ / fitness_);
}

void Gland::unit_test() {
    std::cerr << __PRETTY_FUNCTION__ << std::endl;
    std::cerr.precision(15);
}
