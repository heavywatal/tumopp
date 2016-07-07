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

namespace tumopp {

double Cell::BIRTH_RATE_;
double Cell::DEATH_RATE_;
double Cell::MIGRATION_RATE_;
double Cell::GAMMA_SHAPE_;
double Cell::PROB_SYMMETRIC_DIVISION_;
size_t Cell::MAX_PROLIFERATION_CAPACITY_;
double Cell::MUTATION_RATE_;
double Cell::DRIVER_RATE_BIRTH_;
double Cell::DRIVER_RATE_DEATH_;
double Cell::DRIVER_RATE_MIGRA_;
double Cell::DRIVER_MEAN_BIRTH_;
double Cell::DRIVER_MEAN_DEATH_;
double Cell::DRIVER_MEAN_MIGRA_;
double Cell::DRIVER_SD_BIRTH_;
double Cell::DRIVER_SD_DEATH_;
double Cell::DRIVER_SD_MIGRA_;

namespace {
    std::bernoulli_distribution BERN_BIRTH;
    std::bernoulli_distribution BERN_DEATH;
    std::bernoulli_distribution BERN_MIGRA;
    std::normal_distribution<double> GAUSS_BIRTH;
    std::normal_distribution<double> GAUSS_DEATH;
    std::normal_distribution<double> GAUSS_MIGRA;
}

//! Program options
/*! @return Program options description

    Command line option | Symbol                  | Variable
    --------------------| ----------------------- | -------------------------
    `-b,--beta0`        | \f$\beta_0\f$           | Cell::BIRTH_RATE_
    `-d,--delta0`       | \f$\delta_0\f$          | Cell::DEATH_RATE_
    `-m,--rho0`         | \f$\rho_0\f$            | Cell::MIGRATION_RATE_
    `-k,--shape`        | \f$k\f$                 | Cell::GAMMA_SHAPE_
    `-p,--symmetric`    | \f$p_s\f$               | Cell::PROB_SYMMETRIC_DIVISION_
    `-r,--prolif`       | \f$\omega_\text{max}\f$ | Cell::MAX_PROLIFERATION_CAPACITY_
    `-u,--mutation`     | \f$\mu\f$               | Cell::MUTATION_RATE_
    `--ub`              |                         | Cell::DRIVER_RATE_BIRTH_
    `--ud`              |                         | Cell::DRIVER_RATE_DEATH_
    `--um`              |                         | Cell::DRIVER_RATE_MIGRA_
    `--mb`              |                         | Cell::DRIVER_MEAN_BIRTH_
    `--md`              |                         | Cell::DRIVER_MEAN_DEATH_
    `--mm`              |                         | Cell::DRIVER_MEAN_MIGRA_
    `--sb`              |                         | Cell::DRIVER_SD_BIRTH_
    `--sd`              |                         | Cell::DRIVER_SD_DEATH_
    `--sm`              |                         | Cell::DRIVER_SD_MIGRA_
*/
boost::program_options::options_description Cell::opt_description() {
    namespace po = boost::program_options;
    po::options_description desc{"Cell"};
    desc.add_options()
        ("beta0,b", po::value<double>(&BIRTH_RATE_)->default_value(1.0))
        ("delta0,d", po::value<double>(&DEATH_RATE_)->default_value(0.0))
        ("rho0,m", po::value<double>(&MIGRATION_RATE_)->default_value(0.0))
        ("shape,k", po::value<double>(&GAMMA_SHAPE_)->default_value(1.0))
        ("symmetric,p", po::value<double>(&PROB_SYMMETRIC_DIVISION_)->default_value(1.0))
        ("prolif,r", po::value<size_t>(&MAX_PROLIFERATION_CAPACITY_)->default_value(10))
        ("mutation,u", po::value<double>(&MUTATION_RATE_)->default_value(1e-1))
        ("ub", po::value<double>(&DRIVER_RATE_BIRTH_)->default_value(0.0))
        ("ud", po::value<double>(&DRIVER_RATE_DEATH_)->default_value(0.0))
        ("um", po::value<double>(&DRIVER_RATE_MIGRA_)->default_value(0.0))
        ("mb", po::value<double>(&DRIVER_MEAN_BIRTH_)->default_value(0.0))
        ("md", po::value<double>(&DRIVER_MEAN_DEATH_)->default_value(0.0))
        ("mm", po::value<double>(&DRIVER_MEAN_MIGRA_)->default_value(0.0))
        ("sb", po::value<double>(&DRIVER_SD_BIRTH_)->default_value(0.0))
        ("sd", po::value<double>(&DRIVER_SD_DEATH_)->default_value(0.0))
        ("sm", po::value<double>(&DRIVER_SD_MIGRA_)->default_value(0.0))
    ;
    return desc;
}

void Cell::init_distributions() {
    BERN_BIRTH = std::bernoulli_distribution(DRIVER_RATE_BIRTH_);
    BERN_DEATH = std::bernoulli_distribution(DRIVER_RATE_DEATH_);
    BERN_MIGRA = std::bernoulli_distribution(DRIVER_RATE_MIGRA_);
    GAUSS_BIRTH = std::normal_distribution<double>(DRIVER_MEAN_BIRTH_, DRIVER_SD_BIRTH_);
    GAUSS_DEATH = std::normal_distribution<double>(DRIVER_MEAN_DEATH_, DRIVER_SD_DEATH_);
    GAUSS_MIGRA = std::normal_distribution<double>(DRIVER_MEAN_MIGRA_, DRIVER_SD_MIGRA_);
}

Cell::Cell(const Cell& other):
    coord_(other.coord_),
    birth_rate_(other.birth_rate_),
    death_rate_(other.death_rate_),
    migra_rate_(other.migra_rate_),
    type_(other.type_),
    proliferation_capacity_(other.proliferation_capacity_),
    genealogy_(other.genealogy_) {
    if (type_ == CellType::stem) {
        if (!std::bernoulli_distribution(PROB_SYMMETRIC_DIVISION_)(wtl::sfmt())) {
            type_ = CellType::nonstem;
        }
    }
}

std::string Cell::mutate() {
    auto oss = wtl::make_oss();
    if (BERN_BIRTH(wtl::sfmt())) {
        auto s = GAUSS_BIRTH(wtl::sfmt());
        oss << genealogy_.back() << "\tbirth\t" << s << "\n";
        birth_rate_ *= (s += 1.0);
    }
    if (BERN_DEATH(wtl::sfmt())) {
        auto s = GAUSS_DEATH(wtl::sfmt());
        oss << genealogy_.back() << "\tdeath\t" << s << "\n";
        death_rate_ *= (s += 1.0);
    }
    if (BERN_MIGRA(wtl::sfmt())) {
        auto s = GAUSS_MIGRA(wtl::sfmt());
        oss << genealogy_.back() << "\tmigra\t" << s << "\n";
        migra_rate_ *= (s += 1.0);
    }
    return oss.str();
}

double Cell::delta_time(const double positional_value) {
    double t_birth = std::numeric_limits<double>::infinity();
    if (proliferation_capacity_ > 0) {
        double mu = 1.0;
        mu /= birth_rate_;
        mu /= positional_value;
        mu -= elapsed_;
        double theta = std::max(mu / GAMMA_SHAPE_, 0.0);
        std::gamma_distribution<double> gamma(GAMMA_SHAPE_, theta);
        t_birth = gamma(wtl::sfmt());
    }
    std::exponential_distribution<double> exponential_death(death_rate_);
    std::exponential_distribution<double> exponential_migra(migra_rate_);
    const double t_death = exponential_death(wtl::sfmt());
    const double t_migra = exponential_migra(wtl::sfmt());

    if (t_birth < t_death && t_birth < t_migra) {
        next_event_ = Event::birth;
        elapsed_ = 0.0;
        return t_birth;
    } else if (t_death < t_migra) {
        next_event_ = Event::death;
        return t_death;
    } else {
        next_event_ = Event::migration;
        elapsed_ = t_birth - t_migra;
        return t_migra;
    }
}

std::vector<int> Cell::has_mutations_of(const std::vector<size_t>& mutants) {
    std::vector<int> genotype;
    genotype.reserve(mutants.size());
    for (const size_t mut: mutants) {
        if (std::find(genealogy_.begin(), genealogy_.end(), mut) != genealogy_.end()) {
            genotype.push_back(1);
        } else {
            genotype.push_back(0);
        }
    }
    return genotype;
}

size_t Cell::branch_length(const Cell& other) const {
    const size_t this_len = genealogy_.size();
    const size_t other_len = other.genealogy_.size();
    const size_t shorter = std::min(this_len, other_len);
    size_t branch = this_len + other_len;
    for (size_t i=0; i<shorter; ++i) {
        if (genealogy_[i] != other.genealogy_[i]) break;
        branch -= 2;
    }
    return branch;
}

std::string Cell::header(const char* sep) {
    std::ostringstream oss;
    oss << "x" << sep << "y" << sep << "z" << sep
        << "genealogy" << sep
        << "birth" << sep << "death" << sep
        << "beta" << sep << "delta" << sep << "rho";
    return oss.str();
}

std::ostream& Cell::write(std::ostream& ost, const char* sep) const {
    int z = 0;
    if (coord_.size() > 2) {z = coord_[2];}
    return ost
        << coord_[0] << sep << coord_[1] << sep << z << sep
        << wtl::join(genealogy_, ":") << sep
        << time_of_birth_ << sep << time_of_death_ << sep
        << birth_rate_ << sep
        << death_rate_ << sep
        << migra_rate_;
}

std::string Cell::str(const char* sep) const {
    auto oss = wtl::make_oss();
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
    std::cerr << Cell::header("\t") << "\n" << cell << std::endl;
}

} // namespace tumopp
