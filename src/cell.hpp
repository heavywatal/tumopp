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

namespace boost {namespace program_options {class options_description;}}

namespace tumopp {

//! C1 cell type
enum class CellType: uint_fast8_t {
   stem,
   nonstem,
};

//! event types
enum class Event: uint_fast8_t {
   birth,
   death,
   migration,
};

/*! @brief Parameters for the probability distributions of waiting times
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

/*! @brief Parameters for driver mutations
*/
struct DriverParams {
    //! \f$\mu_\beta\f$
    double RATE_BIRTH = 0.0;
    //! \f$\mu_\delta\f$
    double RATE_DEATH = 0.0;
    //! \f$\mu_\rho\f$
    double RATE_MIGRA = 0.0;
    //! \f$\bar s_\beta\f$
    double MEAN_BIRTH = 0.0;
    //! \f$\bar s_\delta\f$
    double MEAN_DEATH = 0.0;
    //! \f$\bar s_\rho\f$
    double MEAN_MIGRA = 0.0;
    //! \f$\sigma_\beta\f$
    double SD_BIRTH = 0.0;
    //! \f$\sigma_\delta\f$
    double SD_DEATH = 0.0;
    //! \f$\sigma_\rho\f$
    double SD_MIGRA = 0.0;
};

/*! @brief Cancer cell
*/
class Cell {
  public:
    //! Default constructor
    Cell() = default;
    //! Constructor for first cells
    Cell(const std::valarray<int>& v, unsigned i=0,
         std::shared_ptr<EventRates> er=std::make_shared<EventRates>()) noexcept:
      coord_(v), event_rates_(er), id_(i) {}
    //! Copy constructor
    Cell(const Cell& other) noexcept;
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
    double delta_time(double positional_value);

    //! Check mutant ancestors
    std::vector<unsigned> has_mutations_of(const std::vector<unsigned>&) const;

    //! Branch length (# of divisions) between two cells
    size_t branch_length(const Cell&) const;

    //! Set #time_of_birth_; reset other properties
    void set_time_of_birth(double t, unsigned i, const std::shared_ptr<Cell>& ancestor) noexcept {
        time_of_birth_ = t;
        id_ = i;
        ancestor_ = ancestor;
        if (type_ == CellType::nonstem) {--proliferation_capacity_;}
    }
    //! Set #event_rates_->death_prob and #next_event_
    void set_cycle_dependent_death(double death_prob);
    //! Increase #event_rates_->death_rate to birth_rate() for simulating Moran-like situation
    void increase_death_rate() noexcept {event_rates_->death_rate = birth_rate();}

    //! @name Setter functions
    //@{
    void set_coord(const std::valarray<int>& v) noexcept {coord_ = v;}
    void set_time_of_death(double t) noexcept {time_of_death_ = t;}
    void set_elapsed(double t) noexcept {elapsed_ = t;}
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
    //! TSV
    std::string str() const;
    friend std::ostream& operator<< (std::ostream&, const Cell&);

    //! Initialize probability distributions for events
    static void init_distributions(const DriverParams& dp=DriverParams{});

    static boost::program_options::options_description opt_description();

  private:
    //! Accumulate ancestral #id_
    std::unordered_set<unsigned> traceback() const;
    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
    // Data member

    //! \f$k\f$
    static double GAMMA_SHAPE_;
    //! \f$p_s\f$
    static double PROB_SYMMETRIC_DIVISION_;
    //! \f$\omega_\text{max}\f$
    static unsigned MAX_PROLIFERATION_CAPACITY_;

    //! Position in a tumor
    std::valarray<int> coord_;
    //! Set of event rates (copy-on-write)
    std::shared_ptr<EventRates> event_rates_;

    //! C1 cell type
    CellType type_ = CellType::stem;
    //! \f$\omega\f$
    uint_fast8_t proliferation_capacity_ = static_cast<uint_fast8_t>(MAX_PROLIFERATION_CAPACITY_);

    //! next event: birth, death, or migration
    Event next_event_ = Event::birth;
    //! elapsed time by migration
    double elapsed_ = 0.0;

    //! ID
    unsigned id_;
    //! ancestor's ID
    std::shared_ptr<Cell> ancestor_;
    //! time of birth
    double time_of_birth_ = 0.0;
    //! time of death
    double time_of_death_ = 0.0;
};

} // namespace tumopp

#endif /* CELL_HPP_ */
