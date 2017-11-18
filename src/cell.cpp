/*! @file cell.cpp
    @brief Implementation of Cell class
*/
#include "cell.hpp"

#include <sfmt.hpp>
#include <wtl/iostr.hpp>

#include <boost/program_options.hpp>

#include <iostream>
#include <sstream>

namespace tumopp {

double Cell::BIRTH_RATE_;
double Cell::DEATH_RATE_;
double Cell::DEATH_PROB_;
double Cell::MIGRATION_RATE_;
double Cell::GAMMA_SHAPE_;
double Cell::PROB_SYMMETRIC_DIVISION_;
unsigned int Cell::MAX_PROLIFERATION_CAPACITY_;
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
/*! @ingroup params
    @return Program options description

    Command line option | Symbol              | Variable
    ------------------- | ------------------- | -------------------------
    `-b,--beta0`        | \f$\beta_0\f$       | Cell::BIRTH_RATE_
    `-d,--delta0`       | \f$\delta_0\f$      | Cell::DEATH_RATE_
    `-a,--alpha0`       | \f$\alpha_0\f$      | Cell::DEATH_PROB_
    `-m,--rho0`         | \f$\rho_0\f$        | Cell::MIGRATION_RATE_
    `-k,--shape`        | \f$k\f$             | Cell::GAMMA_SHAPE_
    `-p,--symmetric`    | \f$p_s\f$           | Cell::PROB_SYMMETRIC_DIVISION_
    `-r,--prolif`       | \f$\omega_{\max}\f$ | Cell::MAX_PROLIFERATION_CAPACITY_
    `-u,--mutation`     | \f$\mu\f$           | Cell::MUTATION_RATE_
    `--ub`              | \f$\mu_\beta\f$     | Cell::DRIVER_RATE_BIRTH_
    `--ud`              | \f$\mu_\delta\f$    | Cell::DRIVER_RATE_DEATH_
    `--um`              | \f$\mu_\rho\f$      | Cell::DRIVER_RATE_MIGRA_
    `--mb`              | \f$\bar s_\beta\f$  | Cell::DRIVER_MEAN_BIRTH_
    `--md`              | \f$\bar s_\delta\f$ | Cell::DRIVER_MEAN_DEATH_
    `--mm`              | \f$\bar s_\rho\f$   | Cell::DRIVER_MEAN_MIGRA_
    `--sb`              | \f$\sigma_\beta\f$  | Cell::DRIVER_SD_BIRTH_
    `--sd`              | \f$\sigma_\delta\f$ | Cell::DRIVER_SD_DEATH_
    `--sm`              | \f$\sigma_\rho\f$   | Cell::DRIVER_SD_MIGRA_
*/
boost::program_options::options_description Cell::opt_description() {
    namespace po = boost::program_options;
    po::options_description desc{"Cell"};
    desc.add_options()
        ("beta0,b", po::value(&BIRTH_RATE_)->default_value(1.0))
        ("delta0,d", po::value(&DEATH_RATE_)->default_value(0.0))
        ("alpha0,a", po::value(&DEATH_PROB_)->default_value(0.0))
        ("rho0,m", po::value(&MIGRATION_RATE_)->default_value(0.0))
        ("shape,k", po::value(&GAMMA_SHAPE_)->default_value(1.0))
        ("symmetric,p", po::value(&PROB_SYMMETRIC_DIVISION_)->default_value(1.0))
        ("prolif,r", po::value(&MAX_PROLIFERATION_CAPACITY_)->default_value(10))
        ("mutation,u", po::value(&MUTATION_RATE_)->default_value(1e-1))
        ("ub", po::value(&DRIVER_RATE_BIRTH_)->default_value(0.0))
        ("ud", po::value(&DRIVER_RATE_DEATH_)->default_value(0.0))
        ("um", po::value(&DRIVER_RATE_MIGRA_)->default_value(0.0))
        ("mb", po::value(&DRIVER_MEAN_BIRTH_)->default_value(0.0))
        ("md", po::value(&DRIVER_MEAN_DEATH_)->default_value(0.0))
        ("mm", po::value(&DRIVER_MEAN_MIGRA_)->default_value(0.0))
        ("sb", po::value(&DRIVER_SD_BIRTH_)->default_value(0.0))
        ("sd", po::value(&DRIVER_SD_DEATH_)->default_value(0.0))
        ("sm", po::value(&DRIVER_SD_MIGRA_)->default_value(0.0))
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
    death_prob_(other.death_prob_),
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
        double s = GAUSS_BIRTH(wtl::sfmt());
        oss << genealogy_.back() << "\tbirth\t" << s << "\n";
        birth_rate_ *= (s += 1.0);
    }
    if (BERN_DEATH(wtl::sfmt())) {
        double s = GAUSS_DEATH(wtl::sfmt());
        oss << genealogy_.back() << "\tdeath\t" << s << "\n";
        death_rate_ *= (s += 1.0);
        death_prob_ *= (s += 1.0);
    }
    if (BERN_MIGRA(wtl::sfmt())) {
        double s = GAUSS_MIGRA(wtl::sfmt());
        oss << genealogy_.back() << "\tmigra\t" << s << "\n";
        migra_rate_ *= (s += 1.0);
    }
    return oss.str();
}

std::string Cell::force_mutate() {
    birth_rate_ *= (1.0 + DRIVER_MEAN_BIRTH_);
    death_rate_ *= (1.0 + DRIVER_MEAN_DEATH_);
    death_prob_ *= (1.0 + DRIVER_MEAN_DEATH_);
    migra_rate_ *= (1.0 + DRIVER_MEAN_MIGRA_);
    const size_t id = genealogy_.back();
    auto oss = wtl::make_oss();
    oss << id << "\tbirth\t" << DRIVER_MEAN_BIRTH_ << "\n"
        << id << "\tdeath\t" << DRIVER_MEAN_DEATH_ << "\n"
        << id << "\tmigra\t" << DRIVER_MEAN_MIGRA_ << "\n";
    return oss.str();
}

void Cell::increase_death_rate() {
    death_rate_ = birth_rate_;
    elapsed_ = 0.0;
}

double Cell::delta_time(const double positional_value) {
    double t_birth = std::numeric_limits<double>::infinity();
    double t_death = std::numeric_limits<double>::infinity();
    double t_migra = std::numeric_limits<double>::infinity();
    if (proliferation_capacity_ > 0) {
        double mu = 1.0;
        mu /= birth_rate_;
        mu /= positional_value;
        mu -= elapsed_;
        double theta = std::max(mu / GAMMA_SHAPE_, 0.0);
        std::gamma_distribution<double> gamma(GAMMA_SHAPE_, theta);
        t_birth = gamma(wtl::sfmt());
    }
    if (death_rate_ > 0.0) {
        std::exponential_distribution<double> exponential(death_rate_);
        t_death = exponential(wtl::sfmt());
    }
    if (migra_rate_ > 0.0) {
        std::exponential_distribution<double> exponential(migra_rate_);
        t_migra = exponential(wtl::sfmt());
    }

    if (t_birth < t_death && t_birth < t_migra) {
        std::bernoulli_distribution bern_death(death_prob_);
        if (bern_death(wtl::sfmt())) {
            next_event_ = Event::death;
        } else {
            next_event_ = Event::birth;
        }
        elapsed_ = 0.0;
        return t_birth;
    } else if (t_death < t_migra) {
        next_event_ = Event::death;
        return t_death;
    } else {
        next_event_ = Event::migration;
        elapsed_ += t_migra;
        return t_migra;
    }
}

std::vector<unsigned int> Cell::has_mutations_of(const std::vector<size_t>& mutants) {
    std::vector<unsigned int> genotype;
    genotype.reserve(mutants.size());
    for (const auto mut: mutants) {
        if (std::find(genealogy_.begin(), genealogy_.end(), mut) != genealogy_.end()) {
            genotype.push_back(1U);
        } else {
            genotype.push_back(0U);
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

std::string Cell::header() {
    std::ostringstream oss;
    oss << "x\ty\tz\t"
        << "genealogy\t"
        << "birth\tdeath\t"
        << "beta\tdelta\talpha\trho\t"
        << "type\tomega";
    return oss.str();
}

std::ostream& Cell::write(std::ostream& ost) const {
    int z = (coord_.size() > 2U) ? coord_[2] : 0;
    return ost
        << coord_[0] << "\t" << coord_[1] << "\t" << z << "\t"
        << wtl::str_join(genealogy_, ":") << "\t"
        << time_of_birth_ << "\t" << time_of_death_ << "\t"
        << birth_rate_ << "\t"
        << death_rate_ << "\t" << death_prob_ << "\t"
        << migra_rate_ << "\t"
        << static_cast<unsigned int>(type_) << "\t"
        << static_cast<unsigned int>(proliferation_capacity_);
}

std::string Cell::str() const {
    auto oss = wtl::make_oss();
    write(oss);
    return oss.str();
}

//! Stream operator for debug print
std::ostream& operator<< (std::ostream& ost, const Cell& x) {
    return x.write(ost);
}

void Cell::test() {
    std::cerr << __PRETTY_FUNCTION__ << std::endl;
    std::cerr.precision(15);
    Cell cell({1, 2, 3});
    std::cerr << Cell::header() << "\n" << cell << std::endl;
    std::exponential_distribution<double> exponential(0.0);
    std::cerr << "exponential(0): " << exponential(wtl::sfmt()) << std::endl;
}

} // namespace tumopp
