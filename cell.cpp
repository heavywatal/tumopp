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
int Cell::DRIVER_EFFECTS_ = 1;
double Cell::DRIVER_FRACTION_ = 0.0;
double Cell::BIRTH_RATE_ = 1.0;
double Cell::DEATH_RATE_ = 0.0;
double Cell::MIGRATION_RATE_ = 0.0;
double Cell::GAMMA_SHAPE_ = 1.0;
double Cell::PROB_SYMMETRIC_DIVISION_ = 1.0;
size_t Cell::MAX_PROLIFERATION_CAPACITY_ = 10;
size_t Cell::ID_TAIL_ = 0;

std::vector<double> Cell::MUTATION_EFFECTS_;
std::vector<size_t> Cell::MUTANT_IDS_;

//! Program options
/*! @return Program options description

    Command line option | Symbol                  | Variable
    --------------------| ----------------------- | -------------------------
    `-u,--mutation`     | \f$\mu\f$               | Cell::MUTATION_RATE_
    `-s,--sigma`        | \f$\sigma\f$            | Cell::MUTATION_SIGMA_
    `-e,--effect`       |                         | Cell::DRIVER_EFFECTS_
    `-f,--fraction`     |                         | Cell::DRIVER_FRACTION_
    `-b,--birth`        | \f$\beta_0\f$           | Cell::BIRTH_RATE_
    `-d,--death`        | \f$\delta_0\f$          | Cell::DEATH_RATE_
    `-m,--migration`    | \f$m\f$                 | Cell::MIGRATION_RATE_
    `-p,--symmetric`    | \f$p_s\f$               | Cell::PROB_SYMMETRIC_DIVISION_
    `-r,--prolif`       | \f$\rho_\mathrm{max}\f$ | Cell::PROB_SYMMETRIC_DIVISION_
    `-k,--shape`        | \f$k\f$                 | Cell::GAMMA_SHAPE_
*/
boost::program_options::options_description& Cell::opt_description() {
    namespace po = boost::program_options;
    static po::options_description desc{"Cell"};
    desc.add_options()
        ("mutation,u", po::value<double>(&MUTATION_RATE_)->default_value(MUTATION_RATE_))
        ("sigma,s", po::value<double>(&MUTATION_SIGMA_)->default_value(MUTATION_SIGMA_))
        ("effect,e", po::value<int>(&DRIVER_EFFECTS_)->default_value(DRIVER_EFFECTS_))
        ("fraction,f", po::value<double>(&DRIVER_FRACTION_)->default_value(DRIVER_FRACTION_))
        ("birth,b", po::value<double>(&BIRTH_RATE_)->default_value(BIRTH_RATE_))
        ("death,d", po::value<double>(&DEATH_RATE_)->default_value(DEATH_RATE_))
        ("migration,m", po::value<double>(&MIGRATION_RATE_)->default_value(MIGRATION_RATE_))
        ("shape,k", po::value<double>(&GAMMA_SHAPE_)->default_value(GAMMA_SHAPE_))
        ("symmetric,p", po::value<double>(&PROB_SYMMETRIC_DIVISION_)->default_value(PROB_SYMMETRIC_DIVISION_))
        ("prolif,r", po::value<size_t>(&MAX_PROLIFERATION_CAPACITY_)->default_value(MAX_PROLIFERATION_CAPACITY_))
    ;
    return desc;
}

Cell::Cell(const Cell& other):
    coord_(other.coord_), sites_(other.sites_),
    birth_rate_(other.birth_rate_), death_rate_(other.death_rate_),
    type_(other.type_), proliferation_capacity_(other.proliferation_capacity_),
    id_(++ID_TAIL_), mother_(other.id_) {
    if (type_ == CellType::stem) {
        if (!std::bernoulli_distribution(PROB_SYMMETRIC_DIVISION_)(wtl::sfmt())) {
            type_ = CellType::nonstem;
        }
    }
}

void Cell::mutate() {
    static std::bernoulli_distribution bernoulli_driver(DRIVER_FRACTION_);
    static std::normal_distribution<double> normal_sigma(0.0, MUTATION_SIGMA_);
    double effect = 0.0;
    if (DRIVER_FRACTION_ > 0.0) {
        if (bernoulli_driver(wtl::sfmt())) {
            effect = MUTATION_SIGMA_;
        }
    } else {
        effect = normal_sigma(wtl::sfmt());
    }
    sites_.push_back(MUTATION_EFFECTS_.size());
    MUTATION_EFFECTS_.push_back(effect);
    MUTANT_IDS_.push_back(id());
    if (DRIVER_EFFECTS_ & 0b01) {
        birth_rate_ *= (effect += 1.0);
    }
    if (DRIVER_EFFECTS_ & 0b10) {
        death_rate_ *= (1.0 - effect);
    }
}

double Cell::delta_time(const double positional_value) {
    static std::exponential_distribution<double> exponential_migra(MIGRATION_RATE_);
    double t_birth = std::numeric_limits<double>::max();
    if (proliferation_capacity_ > 0) {
        double theta = 1.0;
        theta /= birth_rate();
        theta /= positional_value;
        theta /= GAMMA_SHAPE_;
        std::gamma_distribution<double> gamma(GAMMA_SHAPE_, theta);
        t_birth = gamma(wtl::sfmt());
    }
    std::exponential_distribution<double> exponential(death_rate());
    const double t_death = exponential(wtl::sfmt());
    const double t_migra = exponential_migra(wtl::sfmt());

    if (t_birth < t_death && t_birth < t_migra) {
        next_event_ = Event::birth;
        return t_birth;
    } else if (t_death < t_migra) {
        next_event_ = Event::death;
        return t_death;
    } else {
        next_event_ = Event::migration;
        return t_migra;
    }
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
    oss << "id" << sep << "mother" << sep
        << "birth" << sep << "death" << sep;
    std::vector<std::string> axes{"x", "y", "z"};
    axes.resize(dimensions);
    wtl::ost_join(oss, axes, sep) << sep << "sites" << sep
        << "beta" << sep << "delta\n";
    return oss.str();
}

std::ostream& Cell::write(std::ostream& ost, const std::string& sep) const {
    ost << id_ << sep << mother_ << sep
        << time_of_birth_ << sep << time_of_death_ << sep;
    wtl::ost_join(ost, coord(), sep) << sep;
    wtl::ost_join(ost, sites(), ":") << sep
        << birth_rate_ << sep << death_rate_ << "\n";
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
