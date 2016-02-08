// -*- mode: c++; coding: utf-8 -*-
/*! @file cell.cpp
    @brief Implementation of Cell class
*/
#include "cell.h"

#include <iostream>
#include <sstream>

#include <boost/program_options.hpp>

#include <cxxwtils/prandom.hpp>
#include <cxxwtils/iostr.hpp>

double Cell::MUTATION_RATE_ = 1e-1;
double Cell::MUTATION_SIGMA_ = 0.0;
double Cell::BIRTH_RATE_ = 1.0;
double Cell::DEATH_RATE_ = 0.0;
double Cell::GAMMA_SHAPE_ = 1.0;
size_t Cell::ID_TAIL_ = 0;

std::vector<double> Cell::MUTATION_EFFECTS_;
std::vector<size_t> Cell::MUTANT_IDS_;

//! Program options
/*! @return Program options description

    Command line option | Symbol       | Variable
    --------------------| ------------ | -------------------------
    `-u,--mutation`     | \f$\mu\f$    | Cell::MUTATION_RATE_
    `-s,--sigma`        | \f$\sigma\f$ | Cell::MUTATION_SIGMA_
    `-b,--birth`        | \f$b\f$     | Cell::BIRTH_RATE_
    `-d,--death`        | \f$d\f$     | Cell::DEATH_RATE_
    `-k,--shape`        | \f$k\f$     | Cell::GAMMA_SHAPE_
*/
boost::program_options::options_description& Cell::opt_description() {
    namespace po = boost::program_options;
    static po::options_description desc{"Cell"};
    desc.add_options()
        ("mutation,u", po::value<double>(&MUTATION_RATE_)->default_value(MUTATION_RATE_))
        ("sigma,s", po::value<double>(&MUTATION_SIGMA_)->default_value(MUTATION_SIGMA_))
        ("birth,b", po::value<double>(&BIRTH_RATE_)->default_value(BIRTH_RATE_))
        ("death,d", po::value<double>(&DEATH_RATE_)->default_value(DEATH_RATE_))
        ("shape,k", po::value<double>(&GAMMA_SHAPE_)->default_value(GAMMA_SHAPE_))
    ;
    return desc;
}

void Cell::mutate() {
    sites_.push_back(MUTATION_EFFECTS_.size());
    MUTATION_EFFECTS_.push_back(wtl::prandom().gauss(0.0, MUTATION_SIGMA_));
    MUTANT_IDS_.push_back(id());
    fitness_ += MUTATION_EFFECTS_.back();
}

double Cell::delta_time() const {
    double theta = 1.0;
//    theta /= instantaneous_event_rate();
    theta /= (birth_rate() + death_rate());
    theta /= GAMMA_SHAPE_;
    std::gamma_distribution<double> dist(GAMMA_SHAPE_, theta);
    return dist(wtl::sfmt());
    // return 1.0 / (birth_rate() + death_rate());
}

std::vector<size_t> Cell::haplotype(std::vector<size_t> segsites) const {
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

std::string Cell::header(const size_t dimensions, const std::string& sep) {
    std::ostringstream oss;
    oss << "id" << sep << "mother" << sep << "ancestor" << sep
        << "birth" << sep << "death" << sep;
    std::vector<std::string> axes{"x", "y", "z"};
    axes.resize(dimensions);
    wtl::ost_join(oss, axes, sep) << sep << "sites" << sep << "fitness\n";
    return oss.str();
}

std::ostream& Cell::write(std::ostream& ost, const std::string& sep) const {
    ost << id_ << sep << mother_ << sep << ancestor_ << sep
        << time_of_birth_ << sep << time_of_death_ << sep;
    wtl::ost_join(ost, coord(), sep) << sep;
    wtl::ost_join(ost, sites(), ":") << sep
        << fitness() << "\n";
    return ost;
}

//! Stream operator for debug print
std::ostream& operator<< (std::ostream& ost, const Cell& x) {
    return x.write(ost, "\t");
}

void Cell::unit_test() {
    std::cerr << __PRETTY_FUNCTION__ << std::endl;
    std::cerr.precision(15);
    Cell cell({1, 2, 3});
    std::cerr << Cell::header(3, "\t") << "\n" << cell << std::endl;
}
