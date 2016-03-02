// -*- mode: c++; coding: utf-8 -*-
/*! @file cell.h
    @brief Interface of Cell class
*/
#pragma once
#ifndef CELL_H_
#define CELL_H_

#include <cmath>
#include <vector>
#include <string>

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

namespace boost {
    namespace program_options {
        class options_description;
    }
}


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
    Cell(const std::vector<int>& v): coord_(v), id_(++ID_TAIL_), ancestor_(id_) {}
    //! Copy constructor
    Cell(const Cell& other):
        coord_(other.coord_), sites_(other.sites_),
        birth_rate_(other.birth_rate_), death_rate_(other.death_rate_),
        id_(++ID_TAIL_), mother_(other.id_), ancestor_(other.ancestor_) {}
    //! Destructor
    ~Cell() {--ID_TAIL_;}
    //! Copy assignment operator
    Cell& operator=(const Cell&) = delete;
    //! Move constructor
    Cell(Cell&& other) = default;
    //! Move assignment operator
    Cell& operator=(Cell&&) = default;

    //! Mutate and record
    void mutate();

    //! Setter
    void set_coord(const std::vector<int>& v) {coord_ = v;}
    void set_time_of_birth(const double t) {time_of_birth_ = t;}
    void set_time_of_death(const double t) {time_of_death_ = t;}

    //! Getter
    bool is_dividing() const {return next_event_ == Event::birth;}
    bool is_dying() const {return next_event_ == Event::death;}
    bool is_migrating() const {return next_event_ == Event::migration;}
    double mutation_rate() const {return MUTATION_RATE_;}
    double delta_time();

    const std::vector<int>& coord() const {return coord_;}
    const std::vector<size_t>& sites() const {return sites_;}
    size_t id() const {return id_;}
    size_t mother() const {return mother_;}
    size_t ancestor() const {return ancestor_;}
    double time_of_birth() const {return time_of_birth_;}
    double time_of_death() const {return time_of_death_;}

    //! convert site positions to 01 vector
    std::vector<size_t> haplotype(std::vector<size_t> segsites) const;

    static const std::vector<double>& MUTATION_EFFECTS() {return MUTATION_EFFECTS_;}
    static const std::vector<size_t>& MUTANT_IDS() {return MUTANT_IDS_;}

    static std::string header(const size_t dimensions, const std::string& sep);
    std::ostream& write(std::ostream& ost, const std::string& sep) const;
    friend std::ostream& operator<< (std::ostream&, const Cell&);

    //! Unit test
    static void unit_test();
    static boost::program_options::options_description& opt_description();

  private:
    //! finite per capita rates
    double birth_rate() const {return birth_rate_;}
    double death_rate() const {return death_rate_;}
    double growth_rate() const {return 1.0 + birth_rate_ - death_rate_;}
    //! instantaneous rate for time increment
    double instantaneous_event_rate() const {
        const double lambda = growth_rate();
        return (birth_rate() + death_rate()) * std::log(lambda) / (lambda - 1.0);
    }

    //! per cell division
    static double MUTATION_RATE_;

    static double MUTATION_SIGMA_;

    //! 1: birth, 2: death, 3: both
    static int DRIVER_EFFECTS_;

    static double DRIVER_FRACTION_;

    static double BIRTH_RATE_;

    static double DEATH_RATE_;

    static double MIGRATION_RATE_;
    //! k
    static double GAMMA_SHAPE_;

    static size_t ID_TAIL_;

    //! The history of mutation effects
    static std::vector<double> MUTATION_EFFECTS_;
    static std::vector<size_t> MUTANT_IDS_;

    //! Position in a tumor
    std::vector<int> coord_;
    //! Mutated sites (infinite-site model)
    std::vector<size_t> sites_;

    double birth_rate_ = BIRTH_RATE_;
    double death_rate_ = DEATH_RATE_;

    //! Extra data
    size_t id_ = 0;
    size_t mother_ = 0;
    size_t ancestor_ = 0;
    double time_of_birth_ = 0.0;
    double time_of_death_ = 0.0;
    Event next_event_ = Event::birth;
};

#endif /* CELL_H_ */
