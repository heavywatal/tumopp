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
size_t Gland::ID_TAIL_ = 0;

std::vector<double> Gland::MUTATION_EFFECTS_;
std::vector<size_t> Gland::MUTANT_IDS_;

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
    sites_.push_back(MUTATION_EFFECTS_.size());
    MUTATION_EFFECTS_.push_back(wtl::prandom().gauss(0.0, MUTATION_SIGMA_));
    MUTANT_IDS_.push_back(id());
}

bool Gland::bernoulli_mutation() {
    return wtl::prandom().bernoulli(MUTATION_RATE_ * CELLS_PER_GLAND_);
}

bool Gland::bernoulli_birth() const {
    return wtl::prandom().bernoulli(BIRTH_RATE_);
}

bool Gland::bernoulli_death() const {
    return wtl::prandom().bernoulli(DEATH_RATE_);
}

std::vector<size_t> Gland::haplotype(std::vector<size_t> segsites) const {
    size_t i = 0;
    for (auto& x: segsites) {
        if (i >= sites_.size() || x < sites_[i]) {
            x = 0;
        } else {
            x = 1;
            ++i;
        }
    }
    return segsites;
}

std::string Gland::header(const size_t dimensions, const std::string& sep) {
    std::ostringstream oss;
    oss << "id" << sep << "mother" << sep << "ancestor" << sep;
    std::vector<std::string> axes{"x", "y", "z"};
    axes.resize(dimensions);
    wtl::ost_join(oss, axes, sep) << sep << "sites" << sep << "fitness\n";
    return oss.str();
}

std::ostream& Gland::write(std::ostream& ost, const std::string& sep) const {
    ost << id_ << sep << mother_ << sep << ancestor_ << sep;
    wtl::ost_join(ost, coord(), sep) << sep;
    wtl::ost_join(ost, sites(), ":") << sep
        << fitness() << "\n";
    return ost;
}

//! Stream operator for debug print
std::ostream& operator<< (std::ostream& ost, const Gland& gland) {
    return gland.write(ost, "\t");
}

void Gland::unit_test() {
    std::cerr << __PRETTY_FUNCTION__ << std::endl;
    std::cerr.precision(15);
    Gland gland({1, 2, 3});
    std::cerr << Gland::header(3, "\t") << "\n" << gland << std::endl;
}
