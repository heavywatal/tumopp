// -*- mode: c++; coding: utf-8 -*-
/*! @file cell.cpp
    @brief Implementation of Cell class
*/
#include "cell.hpp"

#include <iostream>
#include <sstream>

#include <boost/program_options.hpp>

#include <cxxwtils/prandom.hpp>
#include <cxxwtils/iostr.hpp>

double Cell::MUTATION_RATE_ = 1e-1;
double Cell::DRIVER_RATE_BIRTH_ = 0.0;
double Cell::DRIVER_RATE_DEATH_ = 0.0;
double Cell::DRIVER_RATE_MIGRA_ = 0.0;
double Cell::DRIVER_MEAN_BIRTH_ = 0.0;
double Cell::DRIVER_MEAN_DEATH_ = 0.0;
double Cell::DRIVER_MEAN_MIGRA_ = 0.0;
double Cell::DRIVER_SD_BIRTH_ = 0.0;
double Cell::DRIVER_SD_DEATH_ = 0.0;
double Cell::DRIVER_SD_MIGRA_ = 0.0;
double Cell::BIRTH_RATE_ = 1.0;
double Cell::DEATH_RATE_ = 0.0;
double Cell::MIGRATION_RATE_ = 0.0;
double Cell::GAMMA_SHAPE_ = 1.0;
double Cell::PROB_SYMMETRIC_DIVISION_ = 1.0;
size_t Cell::MAX_PROLIFERATION_CAPACITY_ = 10;

//! Program options
/*! @return Program options description

    Command line option | Symbol                  | Variable
    --------------------| ----------------------- | -------------------------
    `-u,--mutation`     | \f$\mu\f$               | Cell::MUTATION_RATE_
    `-b,--birth`        | \f$\beta_0\f$           | Cell::BIRTH_RATE_
    `-d,--death`        | \f$\delta_0\f$          | Cell::DEATH_RATE_
    `-m,--migration`    | \f$\rho_0\f$            | Cell::MIGRATION_RATE_
    `-p,--symmetric`    | \f$p_s\f$               | Cell::PROB_SYMMETRIC_DIVISION_
    `-r,--prolif`       | \f$\omega_\text{max}\f$ | Cell::PROB_SYMMETRIC_DIVISION_
    `-k,--shape`        | \f$k\f$                 | Cell::GAMMA_SHAPE_
*/
boost::program_options::options_description& Cell::opt_description() {
    namespace po = boost::program_options;
    static po::options_description desc{"Cell"};
    desc.add_options()
        ("mutation,u", po::value<double>(&MUTATION_RATE_)->default_value(MUTATION_RATE_))
        ("ub", po::value<double>(&DRIVER_RATE_BIRTH_)->default_value(DRIVER_RATE_BIRTH_))
        ("ud", po::value<double>(&DRIVER_RATE_DEATH_)->default_value(DRIVER_RATE_DEATH_))
        ("um", po::value<double>(&DRIVER_RATE_MIGRA_)->default_value(DRIVER_RATE_MIGRA_))
        ("mb", po::value<double>(&DRIVER_MEAN_BIRTH_)->default_value(DRIVER_MEAN_BIRTH_))
        ("md", po::value<double>(&DRIVER_MEAN_DEATH_)->default_value(DRIVER_MEAN_DEATH_))
        ("mm", po::value<double>(&DRIVER_MEAN_MIGRA_)->default_value(DRIVER_MEAN_MIGRA_))
        ("sb", po::value<double>(&DRIVER_SD_BIRTH_)->default_value(DRIVER_SD_BIRTH_))
        ("sd", po::value<double>(&DRIVER_SD_DEATH_)->default_value(DRIVER_SD_DEATH_))
        ("sm", po::value<double>(&DRIVER_SD_MIGRA_)->default_value(DRIVER_SD_MIGRA_))
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
    migra_rate_(other.migra_rate_),
    type_(other.type_), proliferation_capacity_(other.proliferation_capacity_),
    ancestors_(other.ancestors_) {
    ancestors_.push_back(other.id_);
    if (type_ == CellType::stem) {
        if (!std::bernoulli_distribution(PROB_SYMMETRIC_DIVISION_)(wtl::sfmt())) {
            type_ = CellType::nonstem;
        }
    }
}

void Cell::mutate() {
    // TODO
}

double Cell::delta_time(const double positional_value) {
    double t_birth = std::numeric_limits<double>::infinity();
    if (proliferation_capacity_ > 0) {
        double theta = 1.0;
        theta /= birth_rate_;
        theta /= positional_value;
        theta /= GAMMA_SHAPE_;
        std::gamma_distribution<double> gamma(GAMMA_SHAPE_, theta);
        t_birth = gamma(wtl::sfmt());
    }
    std::exponential_distribution<double> exponential_death(death_rate_);
    std::exponential_distribution<double> exponential_migra(migra_rate_);
    const double t_death = exponential_death(wtl::sfmt());
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

std::vector<int> Cell::is_descendant_of(const std::vector<size_t>& mutants) {
    std::vector<int> genotype;
    genotype.reserve(mutants.size());
    for (const size_t mut: mutants) {
        if (std::find(ancestors_.begin(), ancestors_.end(), mut) != ancestors_.end()) {
            genotype.push_back(1);
        } else {
            genotype.push_back(0);
        }
    }
    return genotype;
}

std::string Cell::header(const size_t dimensions, const char* sep) {
    std::vector<std::string> axes{"x", "y", "z"};
    axes.resize(dimensions);
    std::ostringstream oss;
    oss << "id" << sep << "ancestors" << sep
        << "birth" << sep << "death" << sep
        << wtl::join(axes, sep) << sep << "sites" << sep
        << "beta" << sep << "delta" << sep << "rho\n";
    return oss.str();
}

std::ostream& Cell::write(std::ostream& ost, const char* sep) const {
    return ost << id_ << sep << wtl::join(ancestors_, ":") << sep
        << time_of_birth_ << sep << time_of_death_ << sep
        << wtl::join(coord(), sep) << sep
        << wtl::join(sites(), ":") << sep
        << birth_rate_ << sep
        << death_rate_ << sep
        << migra_rate_ << "\n";
}

std::string Cell::str(const char* sep) const {
    std::ostringstream oss;
    write(oss, sep);
    return oss.str();
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
