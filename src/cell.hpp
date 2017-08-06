// -*- mode: c++; coding: utf-8 -*-
/*! @file cell.hpp
    @brief Interface of Cell class
*/
#pragma once
#ifndef CELL_HPP_
#define CELL_HPP_

#include <cstdint>
#include <vector>
#include <valarray>
#include <string>
#include <random>

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

namespace boost {
    namespace program_options {
        class options_description;
    }
}

namespace tumopp {

enum class CellType: uint_fast8_t {
   stem,
   nonstem,
};

enum class Event: uint_fast8_t {
   birth,
   death,
   migration,
};


class Cell {
  public:
    //! Default constructor
    Cell() = default;
    //! Constructor for first cells
    Cell(const std::valarray<int>& v, const size_t i=0): coord_(v) {
        genealogy_.push_back(i);
    }
    //! Copy constructor
    Cell(const Cell& other);
    //! Destructor
    ~Cell() = default;
    //! Copy assignment operator
    Cell& operator=(const Cell&) = delete;
    //! Move constructor
    Cell(Cell&& other) = default;
    //! Move assignment operator
    Cell& operator=(Cell&&) = default;

    //! driver mutation
    std::string mutate();
    std::string force_mutate();

    //! Set death rate to simulate Moran-like model
    void increase_death_rate();

    //! Calc dt and set next_event_
    double delta_time(const double positional_value);

    std::vector<unsigned int> has_mutations_of(const std::vector<size_t>&);

    //! Branch length (# of divisions) between two cells
    size_t branch_length(const Cell&) const;

    //! Setter
    void set_coord(const std::valarray<int>& v) {coord_ = v;}
    void set_time_of_birth(const double t, const size_t i) {
        time_of_birth_ = t;
        genealogy_.push_back(i);
        if (type_ == CellType::nonstem) {--proliferation_capacity_;}
        frustration_ = 0;
    }
    void set_time_of_death(const double t) {time_of_death_ = t;}
    uint_fast8_t frustration() {return ++frustration_;}
    void set_birth_rate(const double x) {birth_rate_ = x;}
    void set_death_rate(const double x) {death_rate_ = x;}
    void set_migra_rate(const double x) {migra_rate_ = x;}

    //! Getter
    double birth_rate() const {return birth_rate_;}
    double death_rate() const {return death_rate_;}
    double migra_rate() const {return migra_rate_;}
    Event next_event() const {return next_event_;}
    const std::valarray<int>& coord() const {return coord_;}
    static double MUTATION_RATE() {return MUTATION_RATE_;}

    //! Print
    static std::string header();
    std::ostream& write(std::ostream& ost) const;
    std::string str() const;
    friend std::ostream& operator<< (std::ostream&, const Cell&);

    //! Initialization
    static boost::program_options::options_description opt_description();
    static void init_distributions();

    //! Unit test
    static void unit_test();

  private:
    //! \f$\beta_0\f$
    static double BIRTH_RATE_;
    //! \f$\delta_0\f$
    static double DEATH_RATE_;
    //! \f$\alpha_0\f$, death probability on cell division attempt
    static double DEATH_PROB_;
    //! \f$\rho_0\f$
    static double MIGRATION_RATE_;
    //! \f$k\f$
    static double GAMMA_SHAPE_;
    //! \f$p_s\f$
    static double PROB_SYMMETRIC_DIVISION_;
    //! \f$\omega_\text{max}\f$
    static unsigned int MAX_PROLIFERATION_CAPACITY_;

    //! \f$\mu\f$ per cell division
    static double MUTATION_RATE_;
    //! \f$\mu_\beta\f$
    static double DRIVER_RATE_BIRTH_;
    //! \f$\mu_\delta\f$
    static double DRIVER_RATE_DEATH_;
    //! \f$\mu_\rho\f$
    static double DRIVER_RATE_MIGRA_;
    //! \f$\bar s_\beta\f$
    static double DRIVER_MEAN_BIRTH_;
    //! \f$\bar s_\delta\f$
    static double DRIVER_MEAN_DEATH_;
    //! \f$\bar s_\rho\f$
    static double DRIVER_MEAN_MIGRA_;
    //! \f$\sigma_\beta\f$
    static double DRIVER_SD_BIRTH_;
    //! \f$\sigma_\delta\f$
    static double DRIVER_SD_DEATH_;
    //! \f$\sigma_\rho\f$
    static double DRIVER_SD_MIGRA_;

    //! Position in a tumor
    std::valarray<int> coord_;

    //! \f$\beta\f$
    double birth_rate_ = BIRTH_RATE_;
    //! \f$\delta\f$
    double death_rate_ = DEATH_RATE_;
    //! \f$\alpha\f$
    double death_prob_ = DEATH_PROB_;
    //! \f$\rho\f$
    double migra_rate_ = MIGRATION_RATE_;

    //! C1 cell type
    CellType type_ = CellType::stem;
    //! \f$\omega\f$
    uint_fast8_t proliferation_capacity_ = MAX_PROLIFERATION_CAPACITY_;

    //! next event: birth, death, or migration
    Event next_event_ = Event::birth;
    //! elapsed time by migration
    double elapsed_ = 0.0;
    //! become quiescent after some failures
    uint_fast8_t frustration_ = 0;

    //! ancestor IDs
    std::vector<size_t> genealogy_;
    //! time of birth
    double time_of_birth_ = 0.0;
    //! time of death
    double time_of_death_ = 0.0;
};

} // namespace tumopp

#endif /* CELL_HPP_ */
