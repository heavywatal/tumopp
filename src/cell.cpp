/*! @file cell.cpp
    @brief Implementation of Cell class
*/
#include "cell.hpp"

#include <wtl/iostr.hpp>
#include <sfmt.hpp>
#include <boost/program_options.hpp>

#include <type_traits>

namespace tumopp {

double Cell::BIRTH_RATE_ = 1.0;
double Cell::DEATH_RATE_ = 0.0;
double Cell::DEATH_PROB_ = 0.0;
double Cell::MIGRATION_RATE_ = 0.0;
double Cell::GAMMA_SHAPE_ = 1.0;
double Cell::PROB_SYMMETRIC_DIVISION_ = 1.0;
unsigned int Cell::MAX_PROLIFERATION_CAPACITY_ = 10u;
double Cell::DRIVER_RATE_BIRTH_ = 0.0;
double Cell::DRIVER_RATE_DEATH_ = 0.0;
double Cell::DRIVER_RATE_MIGRA_ = 0.0;
double Cell::DRIVER_MEAN_BIRTH_ = 0.0;
double Cell::DRIVER_MEAN_DEATH_ = 0.0;
double Cell::DRIVER_MEAN_MIGRA_ = 0.0;
double Cell::DRIVER_SD_BIRTH_ = 0.0;
double Cell::DRIVER_SD_DEATH_ = 0.0;
double Cell::DRIVER_SD_MIGRA_ = 0.0;
double Cell::MUTATION_RATE_ = 1e-1;
bool Cell::HAS_AT_LEAST_1_MUTATION_PER_DIVISION_ = false;

namespace {
    std::bernoulli_distribution BERN_BIRTH(0.0);
    std::bernoulli_distribution BERN_DEATH(0.0);
    std::bernoulli_distribution BERN_MIGRA(0.0);
    std::normal_distribution<double> GAUSS_BIRTH(0.0, 0.0);
    std::normal_distribution<double> GAUSS_DEATH(0.0, 0.0);
    std::normal_distribution<double> GAUSS_MIGRA(0.0, 0.0);
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
    `--ub`              | \f$\mu_\beta\f$     | Cell::DRIVER_RATE_BIRTH_
    `--ud`              | \f$\mu_\delta\f$    | Cell::DRIVER_RATE_DEATH_
    `--um`              | \f$\mu_\rho\f$      | Cell::DRIVER_RATE_MIGRA_
    `--mb`              | \f$\bar s_\beta\f$  | Cell::DRIVER_MEAN_BIRTH_
    `--md`              | \f$\bar s_\delta\f$ | Cell::DRIVER_MEAN_DEATH_
    `--mm`              | \f$\bar s_\rho\f$   | Cell::DRIVER_MEAN_MIGRA_
    `--sb`              | \f$\sigma_\beta\f$  | Cell::DRIVER_SD_BIRTH_
    `--sd`              | \f$\sigma_\delta\f$ | Cell::DRIVER_SD_DEATH_
    `--sm`              | \f$\sigma_\rho\f$   | Cell::DRIVER_SD_MIGRA_
    `-u,--mutation`     | \f$\mu\f$           | Cell::MUTATION_RATE_
    `--ms1mut`          |                     | Cell::HAS_AT_LEAST_1_MUTATION_PER_DIVISION_
*/
boost::program_options::options_description Cell::opt_description() {
    namespace po = boost::program_options;
    po::options_description desc{"Cell"};
    desc.add_options()
        ("beta0,b", po::value(&BIRTH_RATE_)->default_value(BIRTH_RATE_))
        ("delta0,d", po::value(&DEATH_RATE_)->default_value(DEATH_RATE_))
        ("alpha0,a", po::value(&DEATH_PROB_)->default_value(DEATH_PROB_))
        ("rho0,m", po::value(&MIGRATION_RATE_)->default_value(MIGRATION_RATE_))
        ("shape,k", po::value(&GAMMA_SHAPE_)->default_value(GAMMA_SHAPE_))
        ("symmetric,p", po::value(&PROB_SYMMETRIC_DIVISION_)->default_value(PROB_SYMMETRIC_DIVISION_))
        ("prolif,r", po::value(&MAX_PROLIFERATION_CAPACITY_)->default_value(MAX_PROLIFERATION_CAPACITY_))
        ("ub", po::value(&DRIVER_RATE_BIRTH_)->default_value(DRIVER_RATE_BIRTH_))
        ("ud", po::value(&DRIVER_RATE_DEATH_)->default_value(DRIVER_RATE_DEATH_))
        ("um", po::value(&DRIVER_RATE_MIGRA_)->default_value(DRIVER_RATE_MIGRA_))
        ("mb", po::value(&DRIVER_MEAN_BIRTH_)->default_value(DRIVER_MEAN_BIRTH_))
        ("md", po::value(&DRIVER_MEAN_DEATH_)->default_value(DRIVER_MEAN_DEATH_))
        ("mm", po::value(&DRIVER_MEAN_MIGRA_)->default_value(DRIVER_MEAN_MIGRA_))
        ("sb", po::value(&DRIVER_SD_BIRTH_)->default_value(DRIVER_SD_BIRTH_))
        ("sd", po::value(&DRIVER_SD_DEATH_)->default_value(DRIVER_SD_DEATH_))
        ("sm", po::value(&DRIVER_SD_MIGRA_)->default_value(DRIVER_SD_MIGRA_))
        ("mutation,u", po::value(&MUTATION_RATE_)->default_value(MUTATION_RATE_))
        ("ms1mut", po::bool_switch(&HAS_AT_LEAST_1_MUTATION_PER_DIVISION_))
    ;
    return desc;
}

void Cell::init_distributions() {
    BERN_BIRTH.param(decltype(BERN_BIRTH)::param_type(DRIVER_RATE_BIRTH_));
    BERN_DEATH.param(decltype(BERN_DEATH)::param_type(DRIVER_RATE_DEATH_));
    BERN_MIGRA.param(decltype(BERN_MIGRA)::param_type(DRIVER_RATE_MIGRA_));
    GAUSS_BIRTH.param(decltype(GAUSS_BIRTH)::param_type(DRIVER_MEAN_BIRTH_, DRIVER_SD_BIRTH_));
    GAUSS_DEATH.param(decltype(GAUSS_DEATH)::param_type(DRIVER_MEAN_DEATH_, DRIVER_SD_DEATH_));
    GAUSS_MIGRA.param(decltype(GAUSS_MIGRA)::param_type(DRIVER_MEAN_MIGRA_, DRIVER_SD_MIGRA_));
}

static_assert(std::is_nothrow_copy_constructible<Cell>{}, "");
static_assert(std::is_nothrow_move_constructible<Cell>{}, "");

Cell::Cell(const Cell& other) noexcept:
    coord_(other.coord_),
    birth_rate_(other.birth_rate_),
    death_rate_(other.death_rate_),
    death_prob_(other.death_prob_),
    migra_rate_(other.migra_rate_),
    type_(other.type_),
    proliferation_capacity_(other.proliferation_capacity_),
    id_(other.id_),
    ancestor_(other.ancestor_),
    time_of_birth_(other.time_of_birth_) {
    if (type_ == CellType::stem) {
        if (!std::bernoulli_distribution(PROB_SYMMETRIC_DIVISION_)(wtl::sfmt64())) {
            type_ = CellType::nonstem;
        }
    }
}

std::string Cell::mutate() {
    auto oss = wtl::make_oss();
    if (BERN_BIRTH(wtl::sfmt64())) {
        double s = GAUSS_BIRTH(wtl::sfmt64());
        oss << id_ << "\tbirth\t" << s << "\n";
        birth_rate_ *= (s += 1.0);
    }
    if (BERN_DEATH(wtl::sfmt64())) {
        double s = GAUSS_DEATH(wtl::sfmt64());
        oss << id_ << "\tdeath\t" << s << "\n";
        death_rate_ *= (s += 1.0);
        death_prob_ *= (s += 1.0);
    }
    if (BERN_MIGRA(wtl::sfmt64())) {
        double s = GAUSS_MIGRA(wtl::sfmt64());
        oss << id_ << "\tmigra\t" << s << "\n";
        migra_rate_ *= (s += 1.0);
    }
    return oss.str();
}

std::string Cell::force_mutate() {
    birth_rate_ *= (1.0 + DRIVER_MEAN_BIRTH_);
    death_rate_ *= (1.0 + DRIVER_MEAN_DEATH_);
    death_prob_ *= (1.0 + DRIVER_MEAN_DEATH_);
    migra_rate_ *= (1.0 + DRIVER_MEAN_MIGRA_);
    auto oss = wtl::make_oss();
    oss << id_ << "\tbirth\t" << DRIVER_MEAN_BIRTH_ << "\n"
        << id_ << "\tdeath\t" << DRIVER_MEAN_DEATH_ << "\n"
        << id_ << "\tmigra\t" << DRIVER_MEAN_MIGRA_ << "\n";
    return oss.str();
}

void Cell::increase_death_rate() noexcept {
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
        t_birth = gamma(wtl::sfmt64());
    }
    if (death_rate_ > 0.0) {
        std::exponential_distribution<double> exponential(death_rate_);
        t_death = exponential(wtl::sfmt64());
    }
    if (migra_rate_ > 0.0) {
        std::exponential_distribution<double> exponential(migra_rate_);
        t_migra = exponential(wtl::sfmt64());
    }

    if (t_birth < t_death && t_birth < t_migra) {
        std::bernoulli_distribution bern_death(death_prob_);
        if (bern_death(wtl::sfmt64())) {
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

std::unordered_set<uint_fast32_t> Cell::traceback() const {
    std::unordered_set<uint_fast32_t> genealogy;
    genealogy.emplace(id_);
    for (std::shared_ptr<Cell> p = ancestor_; p; p = p->ancestor_) {
        genealogy.emplace(p->id_);
    }
    return genealogy;
}

std::vector<uint_fast32_t> Cell::has_mutations_of(const std::vector<uint_fast32_t>& mutants) const {
    const auto genealogy = traceback();
    std::vector<uint_fast32_t> genotype;
    genotype.reserve(mutants.size());
    for (const auto mut: mutants) {
        if (genealogy.find(mut) != genealogy.end()) {
            genotype.push_back(1u);
        } else {
            genotype.push_back(0u);
        }
    }
    return genotype;
}

size_t Cell::branch_length(const Cell& other) const {
    if (id_ == other.id_) return 0u;
    size_t length = 2u;
    uint_fast32_t mrca = 1u;
    const auto genealogy = traceback();
    for (std::shared_ptr<Cell> p = other.ancestor_; p; p = p->ancestor_) {
        if (genealogy.find(p->id_) != genealogy.end()) {
            mrca = p->id_;
            break;
        }
        ++length;
    }
    for (std::shared_ptr<Cell> p = ancestor_; p->id_ > mrca; p = p->ancestor_) {
        ++length;
    }
    return length;
}

std::string Cell::header() {
    std::ostringstream oss;
    oss << "x\ty\tz\t"
        << "id\tancestor\t"
        << "birth\tdeath\t"
        << "beta\tdelta\talpha\trho\t"
        << "type\tomega";
    return oss.str();
}

std::ostream& Cell::write(std::ostream& ost) const {
    int z = (coord_.size() > 2U) ? coord_[2] : 0;
    return ost
        << coord_[0] << "\t" << coord_[1] << "\t" << z << "\t"
        << id_ << "\t"
        << (ancestor_ ? ancestor_->id_ : 0u) << "\t"
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

} // namespace tumopp
