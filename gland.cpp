// -*- mode: c++; coding: utf-8 -*-
/*! @file gland.cpp
    @brief Implementation of Gland class
*/
#include "gland.h"

#include <iostream>
#include <sstream>

#include <boost/program_options.hpp>

#include <cxxwtils/prandom.hpp>
#include <cxxwtils/iostr.hpp>

double Gland::CELLS_PER_GLAND_ = 1e4;
double Gland::MUTATION_RATE_ = 1e-5;
double Gland::MUTATION_SIGMA_ = 0.0;
double Gland::BIRTH_RATE_ = 1.0;
double Gland::DEATH_RATE_ = 0.0;

std::vector<double> Gland::MUTATION_EFFECTS_;

//! Program options
/*! @return Program options description

    Command line option | Symbol       | Variable
    --------------------| ------------ | -------------------------
    `-n,--cells`        | -            | Gland::CELLS_PER_GLAND_
    `-u,--mutation`     | \f$\mu\f$    | Gland::MUTATION_RATE_
    `-s,--sigma`        | \f$\sigma\f$ | Gland::MUTATION_SIGMA_
    `-b,--birth`        | \f$\b\f$     | Gland::BIRTH_RATE_
    `-d,--death`        | \f$\d\f$     | Gland::DEATH_RATE_
*/
boost::program_options::options_description& Gland::opt_description() {
    namespace po = boost::program_options;
    static po::options_description desc{"Gland"};
    desc.add_options()
        ("cells,n", po::value<double>(&CELLS_PER_GLAND_)->default_value(CELLS_PER_GLAND_))
        ("mutation,u", po::value<double>(&MUTATION_RATE_)->default_value(MUTATION_RATE_))
        ("sigma,s", po::value<double>(&MUTATION_SIGMA_)->default_value(MUTATION_SIGMA_))
        ("birth,b", po::value<double>(&BIRTH_RATE_)->default_value(BIRTH_RATE_))
        ("death,d", po::value<double>(&DEATH_RATE_)->default_value(DEATH_RATE_))
    ;
    return desc;
}

void Gland::mutate() {
    MUTATION_EFFECTS_.push_back(wtl::prandom().gauss(0.0, MUTATION_SIGMA_));
    sites_.push_back(MUTATION_EFFECTS_.size());
}

bool Gland::bernoulli_mutation() {
    return wtl::prandom().bernoulli(MUTATION_RATE_ * CELLS_PER_GLAND_);
}

bool Gland::bernoulli_birth() const {
    return wtl::prandom().bernoulli(BIRTH_RATE_);
}

bool Gland::bernoulli_death() const {
    return wtl::prandom().bernoulli(DEATH_RATE_ / fitness());
}

//! Stream operator for debug print
std::ostream& operator<< (std::ostream& ost, const Gland& gland) {
    return ost << gland.sites_;
}

void Gland::unit_test() {
    std::cerr << __PRETTY_FUNCTION__ << std::endl;
    std::cerr.precision(15);
    Gland gland;
    std::cerr << gland << std::endl;
}
