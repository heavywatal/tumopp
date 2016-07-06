// -*- mode: c++; coding: utf-8 -*-
/*! @file cell.h
    @brief Interface of Cell class
*/
#pragma once
#ifndef CELL_HPP_
#define CELL_HPP_

#include <cmath>
#include <vector>
#include <valarray>
#include <string>

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

namespace boost {
    namespace program_options {
        class options_description;
    }
}

namespace tumopp {

enum class CellType: int {
   stem,
   nonstem,
};

enum class Event: int {
   death     = -1,
   migration =  0,
   birth     =  1,
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

    //! Calc dt and set next_event_
    double delta_time(const double positional_value);

    std::vector<int> has_mutations_of(const std::vector<size_t>&);

    //! Branch length (# of divisions) between two cells
    size_t operator-(const Cell&) const;

    //! Setter
    void set_coord(const std::valarray<int>& v) {coord_ = v;}
    void set_time_of_birth(const double t, const size_t i) {
        time_of_birth_ = t;
        genealogy_.push_back(i);
        if (type_ == CellType::nonstem) {--proliferation_capacity_;}
    }
    void set_time_of_death(const double t) {time_of_death_ = t;}

    //! Getter
    Event next_event() const {return next_event_;}
    const std::valarray<int>& coord() const {return coord_;}
    static double MUTATION_RATE() {return MUTATION_RATE_;}

    //! Print
    static std::string header(const char* sep);
    std::ostream& write(std::ostream& ost, const char* sep) const;
    std::string str(const char* sep) const;
    friend std::ostream& operator<< (std::ostream&, const Cell&);

    //! Unit test
    static void unit_test();
    static boost::program_options::options_description opt_description();

  private:
    //! \f$\beta_0\f$
    static double BIRTH_RATE_;
    //! \f$\delta_0\f$
    static double DEATH_RATE_;
    //! \f$\rho_0\f$
    static double MIGRATION_RATE_;
    //! \f$k\f$
    static double GAMMA_SHAPE_;
    //! \f$p_s\f$
    static double PROB_SYMMETRIC_DIVISION_;
    //! \f$\omega_\text{max}\f$
    static size_t MAX_PROLIFERATION_CAPACITY_;

    //! mutation rate per cell division
    static double MUTATION_RATE_;
    static double DRIVER_RATE_BIRTH_;
    static double DRIVER_RATE_DEATH_;
    static double DRIVER_RATE_MIGRA_;
    static double DRIVER_MEAN_BIRTH_;
    static double DRIVER_MEAN_DEATH_;
    static double DRIVER_MEAN_MIGRA_;
    static double DRIVER_SD_BIRTH_;
    static double DRIVER_SD_DEATH_;
    static double DRIVER_SD_MIGRA_;

    //! Position in a tumor
    std::valarray<int> coord_;

    //! between birth and migration
    double elapsed_ = 0.0;

    double birth_rate_ = BIRTH_RATE_;
    double death_rate_ = DEATH_RATE_;
    double migra_rate_ = MIGRATION_RATE_;

    CellType type_ = CellType::stem;
    size_t proliferation_capacity_ = MAX_PROLIFERATION_CAPACITY_;

    Event next_event_ = Event::birth;

    //! Extra data
    std::vector<size_t> genealogy_;
    double time_of_birth_ = 0.0;
    double time_of_death_ = 0.0;
};

} // namespace tumopp

#endif /* CELL_HPP_ */
