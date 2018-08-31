/*! @file cell.hpp
    @brief Interface of Cell class
*/
#pragma once
#ifndef CELL_HPP_
#define CELL_HPP_

#include <cstdint>
#include <vector>
#include <valarray>
#include <unordered_set>
#include <string>
#include <memory>

namespace tumopp {

//! event types
enum class Event: uint_fast8_t {
   birth,
   death,
   migration,
};

//! @brief Parameters for the probability distributions of waiting times
/*! @ingroup params
*/
struct EventRates {
    //! \f$\beta\f$
    double birth_rate = 1.0;
    //! \f$\delta\f$
    double death_rate = 0.0;
    //! \f$\alpha\f$, death probability on cell division attempt
    double death_prob = 0.0;
    //! \f$\rho\f$
    double migra_rate = 0.0;
};

//! @brief Parameters for Cell class
/*! @ingroup params
*/
struct CellParams {
    //! \f$k\f$
    double GAMMA_SHAPE = 1.0;
    //! \f$p_s\f$
    double PROB_SYMMETRIC_DIVISION = 1.0;
    //! \f$\omega_\text{max}\f$
    unsigned MAX_PROLIFERATION_CAPACITY = 10u;
    //! \f$\mu_\beta\f$
    double RATE_BIRTH = 0.0;
    //! \f$\mu_\delta\f$
    double RATE_DEATH = 0.0;
    //! \f$\mu_\alpha\f$
    double RATE_ALPHA = 0.0;
    //! \f$\mu_\rho\f$
    double RATE_MIGRA = 0.0;
    //! \f$\bar s_\beta\f$
    double MEAN_BIRTH = 0.0;
    //! \f$\bar s_\delta\f$
    double MEAN_DEATH = 0.0;
    //! \f$\bar s_\alpha\f$
    double MEAN_ALPHA = 0.0;
    //! \f$\bar s_\rho\f$
    double MEAN_MIGRA = 0.0;
    //! \f$\sigma_\beta\f$
    double SD_BIRTH = 0.0;
    //! \f$\sigma_\delta\f$
    double SD_DEATH = 0.0;
    //! \f$\sigma_\alpha\f$
    double SD_ALPHA = 0.0;
    //! \f$\sigma_\rho\f$
    double SD_MIGRA = 0.0;
};

/*! @brief Cancer cell
*/
class Cell {
  public:
    //! Alias
    using param_type = CellParams;
    //! Default constructor
    Cell() = default;
    //! Constructor for first cells
    Cell(const std::valarray<int>& v, unsigned i=0,
         std::shared_ptr<EventRates> er=std::make_shared<EventRates>()) noexcept:
      coord_(v), event_rates_(er),
      id_(i) {}
    //! Copy constructor
    Cell(const Cell& other) noexcept:
      coord_(other.coord_),
      event_rates_(other.event_rates_),
      proliferation_capacity_(other.proliferation_capacity_),
      id_(other.id_),
      ancestor_(other.ancestor_),
      time_of_birth_(other.time_of_birth_) {}
    //! Destructor
    ~Cell() noexcept = default;
    //! Copy assignment operator
    Cell& operator=(const Cell&) = delete;
    //! Move constructor
    Cell(Cell&& other) = default;
    //! Move assignment operator
    Cell& operator=(Cell&&) = default;

    //! driver mutation
    std::string mutate();
    //! driver mutation on all traits
    std::string force_mutate();

    //! Calc dt and set #next_event_
    double delta_time(double now, double positional_value, bool surrounded=false);

    //! Check mutant ancestors
    std::vector<unsigned> has_mutations_of(const std::vector<unsigned>&) const;

    //! Branch length (# of divisions) between two cells
    size_t branch_length(const Cell&) const;

    //! Change #proliferation_capacity_ stochastically
    void differentiate();
    //! Set #time_of_birth_; reset other properties
    void set_time_of_birth(double t, unsigned i, const std::shared_ptr<Cell>& ancestor) noexcept {
        time_of_birth_ = t;
        id_ = i;
        ancestor_ = ancestor;
        if (is_differentiated()) {--proliferation_capacity_;}
    }
    //! Set #event_rates_->death_prob and #next_event_
    void set_cycle_dependent_death(double death_prob);
    //! Increase #event_rates_->death_rate to birth_rate() for simulating Moran-like situation
    void increase_death_rate() noexcept {event_rates_->death_rate = birth_rate();}
    //! Check #proliferation_capacity_
    bool is_differentiated() const noexcept {return proliferation_capacity_ >= 0;}

    //! @name Setter functions
    //@{
    void set_coord(const std::valarray<int>& v) noexcept {coord_ = v;}
    void set_time_of_death(double t) noexcept {time_of_death_ = t;}
    //@}

    //! @name Getter functions
    //@{
    double birth_rate() const noexcept {return event_rates_->birth_rate;}
    double death_rate() const noexcept {return event_rates_->death_rate;}
    double death_prob() const noexcept {return event_rates_->death_prob;}
    double migra_rate() const noexcept {return event_rates_->migra_rate;}
    Event next_event() const noexcept {return next_event_;}
    const std::valarray<int>& coord() const noexcept {return coord_;}
    //@}

    //! TSV header
    static std::string header();
    //! TSV
    std::ostream& write(std::ostream& ost) const;
    friend std::ostream& operator<< (std::ostream&, const Cell&);

    //! Set #PARAM_
    static void param(const param_type& p);
    //! Get #PARAM_
    static const param_type& param() {return PARAM_;}

  private:
    //! Accumulate ancestral #id_
    std::unordered_set<unsigned> traceback() const;

    //! Parameters shared among instances
    static param_type PARAM_;

    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
    // Data member

    //! Position in a tumor
    std::valarray<int> coord_;
    //! Set of event rates (copy-on-write)
    std::shared_ptr<EventRates> event_rates_;

    //! \f$\omega\f$; stem cell if negative
    int8_t proliferation_capacity_ = -1;

    //! next event: birth, death, or migration
    Event next_event_ = Event::birth;
    //! ID
    unsigned id_;
    //! pointer to the ancestor
    std::shared_ptr<Cell> ancestor_;
    //! time of birth
    double time_of_birth_ = 0.0;
    //! time of death
    double time_of_death_ = 0.0;
};

} // namespace tumopp

#endif /* CELL_HPP_ */
