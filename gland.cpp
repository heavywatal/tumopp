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

double Gland::CELLS_PER_GLAND_ = 1e4;
double Gland::MUTATION_RATE_ = 1e-5;
double Gland::MUTATION_SIGMA_ = 0.0;
double Gland::APOPTOSIS_RATE_ = 0.2;

size_t Gland::MUTATED_SITES_ = 0;

boost::program_options::options_description& Gland::opt_description() {
    namespace po = boost::program_options;
    static po::options_description desc{"Gland"};
    desc.add_options()
        ("cells,c", po::value<double>(&CELLS_PER_GLAND_)->default_value(CELLS_PER_GLAND_))
        ("mutation,u", po::value<double>(&MUTATION_RATE_)->default_value(MUTATION_RATE_))
        ("sigma,s", po::value<double>(&MUTATION_SIGMA_)->default_value(MUTATION_SIGMA_))
        ("apoptosis,a", po::value<double>(&APOPTOSIS_RATE_)->default_value(APOPTOSIS_RATE_))
    ;
    return desc;
}

void Gland::mutate() {
    sites_.push_back(MUTATED_SITES_++);
    fitness_ += wtl::prandom().gauss(0.0, MUTATION_SIGMA_);
    fitness_ = std::max(fitness_, 1.0);
    fitness_ = std::min(fitness_, 5.0);
}

bool Gland::bernoulli_mutation() const {
    return wtl::prandom().bernoulli(MUTATION_RATE_ * CELLS_PER_GLAND_);
}

bool Gland::bernoulli_apoptosis() const {
    return wtl::prandom().bernoulli(APOPTOSIS_RATE_ / fitness_);
}

void Gland::unit_test() {
    std::cerr << __PRETTY_FUNCTION__ << std::endl;
    std::cerr.precision(15);
}
