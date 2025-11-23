/*! @file cell.cpp
    @brief Implementation of Cell class
*/
#include "cell.hpp"

#include <fmt/format.h>
#include <wtl/random.hpp>

namespace tumopp {

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
namespace {

class GammaFactory {
  public:
    GammaFactory(double k) noexcept: shape_(k) {}
    std::gamma_distribution<double> operator()(double mu) {
        const double theta = std::max(mu / shape_, 0.0);
        return std::gamma_distribution<double>(shape_, theta);
    }
    void param(double k) {shape_ = k;}
  private:
    double shape_;
};

template <class URBG>
inline bool bernoulli(double p, URBG& engine) {
    // consume less URBG when p is set to 0 or 1.
    return p >= 1.0 || (p > 0.0 && wtl::generate_canonical(engine) < p);
}

class bernoulli_distribution {
  public:
    bernoulli_distribution(double p) noexcept: p_(p) {}
    template <class URBG>
    bool operator()(URBG& engine) const {
        return p_ >= 1.0 || (p_ > 0.0 && wtl::generate_canonical(engine) < p_);
    }
    void param(double p) {p_ = p;}
  private:
    double p_;
};

GammaFactory GAMMA_FACTORY(Cell::param().GAMMA_SHAPE);
bernoulli_distribution BERN_SYMMETRIC(Cell::param().PROB_SYMMETRIC_DIVISION);
bernoulli_distribution BERN_MUT_BIRTH(Cell::param().RATE_BIRTH);
bernoulli_distribution BERN_MUT_DEATH(Cell::param().RATE_DEATH);
bernoulli_distribution BERN_MUT_ALPHA(Cell::param().RATE_ALPHA);
bernoulli_distribution BERN_MUT_MIG(Cell::param().RATE_MIG);
std::normal_distribution<double> GAUSS_BIRTH(Cell::param().MEAN_BIRTH, Cell::param().SD_BIRTH);
std::normal_distribution<double> GAUSS_DEATH(Cell::param().MEAN_DEATH, Cell::param().SD_DEATH);
std::normal_distribution<double> GAUSS_ALPHA(Cell::param().MEAN_ALPHA, Cell::param().SD_ALPHA);
std::normal_distribution<double> GAUSS_MIG(Cell::param().MEAN_MIG, Cell::param().SD_MIG);

}// namespace
/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

void Cell::param(const param_type& p) {
    PARAM_ = p;
    GAMMA_FACTORY.param(PARAM_.GAMMA_SHAPE);
    BERN_SYMMETRIC.param(PARAM_.PROB_SYMMETRIC_DIVISION);
    BERN_MUT_BIRTH.param(PARAM_.RATE_BIRTH);
    BERN_MUT_DEATH.param(PARAM_.RATE_DEATH);
    BERN_MUT_ALPHA.param(PARAM_.RATE_ALPHA);
    BERN_MUT_MIG.param(PARAM_.RATE_MIG);
    GAUSS_BIRTH.param(decltype(GAUSS_BIRTH)::param_type(PARAM_.MEAN_BIRTH, PARAM_.SD_BIRTH));
    GAUSS_DEATH.param(decltype(GAUSS_DEATH)::param_type(PARAM_.MEAN_DEATH, PARAM_.SD_DEATH));
    GAUSS_ALPHA.param(decltype(GAUSS_ALPHA)::param_type(PARAM_.MEAN_ALPHA, PARAM_.SD_ALPHA));
    GAUSS_MIG.param(decltype(GAUSS_MIG)::param_type(PARAM_.MEAN_MIG, PARAM_.SD_MIG));
}

void Cell::differentiate(urbg_t& engine) {
    if (is_differentiated()) return;
    if (BERN_SYMMETRIC(engine)) return;
    proliferation_capacity_ = static_cast<int8_t>(PARAM_.MAX_PROLIFERATION_CAPACITY);
}

void Cell::mutate(std::string& buffer, urbg_t& engine) {
    if (BERN_MUT_BIRTH(engine)) {
        event_rates_ = std::make_shared<EventRates>(*event_rates_);
        double s = GAUSS_BIRTH(engine);
        fmt::format_to(std::back_inserter(buffer), "{}\tbeta\t{:.9g}\n", id_, s);
        event_rates_->birth_rate *= (s += 1.0);
    }
    if (BERN_MUT_DEATH(engine)) {
        event_rates_ = std::make_shared<EventRates>(*event_rates_);
        double s = GAUSS_DEATH(engine);
        fmt::format_to(std::back_inserter(buffer), "{}\tdelta\t{:.9g}\n", id_, s);
        event_rates_->death_rate *= (s += 1.0);
    }
    if (BERN_MUT_ALPHA(engine)) {
        event_rates_ = std::make_shared<EventRates>(*event_rates_);
        double s = GAUSS_ALPHA(engine);
        fmt::format_to(std::back_inserter(buffer), "{}\talpha\t{:.9g}\n", id_, s);
        event_rates_->death_prob *= (s += 1.0);
    }
    if (BERN_MUT_MIG(engine)) {
        event_rates_ = std::make_shared<EventRates>(*event_rates_);
        double s = GAUSS_MIG(engine);
        fmt::format_to(std::back_inserter(buffer), "{}\trho\t{:.9g}\n", id_, s);
        event_rates_->migration_rate *= (s += 1.0);
    }
}

void Cell::force_mutate(std::string& buffer, urbg_t& engine) {
    event_rates_ = std::make_shared<EventRates>(*event_rates_);
    const double s_birth = GAUSS_BIRTH(engine);
    const double s_death = GAUSS_DEATH(engine);
    const double s_alpha = GAUSS_ALPHA(engine);
    const double s_migration = GAUSS_MIG(engine);
    event_rates_->birth_rate *= (1.0 + s_birth);
    event_rates_->death_rate *= (1.0 + s_death);
    event_rates_->death_prob *= (1.0 + s_alpha);
    event_rates_->migration_rate *= (1.0 + s_migration);
    if (s_birth != 0.0) {
        fmt::format_to(std::back_inserter(buffer), "{}\tbeta\t{:.9g}\n", id_, s_birth);
    }
    if (s_death != 0.0) {
        fmt::format_to(std::back_inserter(buffer), "{}\tdelta\t{:.9g}\n", id_, s_death);
    }
    if (s_alpha != 0.0) {
        fmt::format_to(std::back_inserter(buffer), "{}\talpha\t{:.9g}\n", id_, s_alpha);
    }
    if (s_migration != 0.0) {
        fmt::format_to(std::back_inserter(buffer), "{}\trho\t{:.9g}\n", id_, s_migration);
    }
}

double Cell::delta_time(urbg_t& engine, const double now, const double positional_value, const bool surrounded) {
    double t_birth = std::numeric_limits<double>::infinity();
    double t_death = std::numeric_limits<double>::infinity();
    double t_migration = std::numeric_limits<double>::infinity();
    if (proliferation_capacity_ != 0) {
        double mu = 1.0;
        mu /= birth_rate();
        mu /= positional_value;
        if (!surrounded) mu -= (now - time_of_birth_);
        t_birth = GAMMA_FACTORY(mu)(engine);
    }
    if (death_rate() > 0.0) {
        std::exponential_distribution<double> exponential(death_rate());
        t_death = exponential(engine);
    }
    if (migration_rate() > 0.0) {
        std::exponential_distribution<double> exponential(migration_rate());
        t_migration = exponential(engine);
    }

    if (t_birth < t_death && t_birth < t_migration) {
        next_event_ = bernoulli(death_prob(), engine)
                      ? Event::death : Event::birth;
        return t_birth;
    } else if (t_death < t_migration) {
        next_event_ = Event::death;
        return t_death;
    } else {
        next_event_ = Event::migration;
        return t_migration;
    }
}

void Cell::set_cycle_dependent_death(urbg_t& engine, const double p) {
    //TODO: reduce redundant copy for susceptible cells
    event_rates_ = std::make_shared<EventRates>(*event_rates_);
    event_rates_->death_prob = p;
    next_event_ = bernoulli(p, engine) ? Event::death : Event::birth;
}

const char* Cell::header() {
    return "x\ty\tz\tid\tancestor\t"
           "birth\tdeath\tomega";
}

void Cell::format_to_back(std::string& buffer) const {
    fmt::format_to(std::back_inserter(buffer),
      "{}\t{}\t{}\t{}\t{}\t{:.9g}\t{:.9g}\t{}\n",
      coord_[0], coord_[1], coord_[2],
      id_,
      (ancestor_ ? ancestor_->id_ : 0),
      time_of_birth_, time_of_death_,
      static_cast<int>(proliferation_capacity_));
}

void Cell::traceback(std::string& buffer, std::unordered_set<int>& done) const {
    format_to_back(buffer);
    if (ancestor_ && done.insert(ancestor_->id_).second) {
        ancestor_->traceback(buffer, done);
    }
}

} // namespace tumopp
