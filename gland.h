// -*- mode: c++; coding: utf-8 -*-
/*! @file gland.h
    @brief Interface of Gland class
*/
#pragma once
#ifndef GLAND_H_
#define GLAND_H_

#include <cmath>
#include <vector>
#include <string>

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

namespace boost {
    namespace program_options {
        class options_description;
    }
}


class Gland {
  public:
    //! Default constructor
    Gland() = delete;
    //! Constructor for first cells
    Gland(const std::vector<int>& v): coord_(v), id_(++ID_TAIL_), ancestor_(id_) {}
    //! Copy constructor
    Gland(const Gland& other):
        coord_(other.coord_), sites_(other.sites_), fitness_(other.fitness_),
        id_(++ID_TAIL_), mother_(other.id_), ancestor_(other.ancestor_) {}
    //! Copy assignment operator
    Gland& operator=(const Gland&) = delete;
    //! Move constructor
    Gland(Gland&& other) = default;
    //! Move assignment operator
    Gland& operator=(Gland&&) = default;

    //! Mutate and record
    void mutate();

    //! Setter
    void set_coord(const std::vector<int>& v) {coord_ = v;}
    void set_time_of_birth(const double t) {time_of_birth_ = t;}
    void set_time_of_death(const double t) {time_of_death_ = t;}

    //! Getter
    //! finite per capita rates
    double birth_rate() const {return BIRTH_RATE_ * fitness_;}
    double death_rate() const {return DEATH_RATE_;}
    double increase_rate() const {return 1.0 + birth_rate() - death_rate();}
    double mutation_rate() const {return MUTATION_RATE_ * CELLS_PER_GLAND_;}
    //! instantaneous rate for time increment
    double instantaneous_event_rate() const {
        const double lambda = increase_rate();
        return (birth_rate() + death_rate()) * std::log(lambda) / (lambda - 1.0);
    }

    const std::vector<int>& coord() const {return coord_;}
    const std::vector<size_t>& sites() const {return sites_;}
    size_t id() const {return id_;}
    size_t mother() const {return mother_;}
    size_t ancestor() const {return ancestor_;}
    double time_of_birth() const {return time_of_birth_;}
    double time_of_death() const {return time_of_death_;}

    //! convert site positions to 01 vector
    std::vector<size_t> haplotype(std::vector<size_t> segsites) const;

    double fitness() const {return fitness_;}

    static const std::vector<double>& MUTATION_EFFECTS() {return MUTATION_EFFECTS_;}
    static const std::vector<size_t>& MUTANT_IDS() {return MUTANT_IDS_;}

    static std::string header(const size_t dimensions, const std::string& sep);
    std::ostream& write(std::ostream& ost, const std::string& sep) const;
    friend std::ostream& operator<< (std::ostream&, const Gland&);

    //! Unit test
    static void unit_test();
    static boost::program_options::options_description& opt_description();

  private:
    //! 80 billion cells / 8 million glands = 10000
    static double CELLS_PER_GLAND_;

    //! per gland per division {1e-8 to 1e-4}
    static double MUTATION_RATE_;

    //! {0, 0.2, 0.6}
    static double MUTATION_SIGMA_;

    static double BIRTH_RATE_;

    static double DEATH_RATE_;

    static size_t ID_TAIL_;

    //! The history of mutation effects
    static std::vector<double> MUTATION_EFFECTS_;
    static std::vector<size_t> MUTANT_IDS_;

    //! Position in a tumor
    std::vector<int> coord_;
    //! Mutated sites (infinite-site model)
    std::vector<size_t> sites_;
    //! Sum of mutation effects + 1
    double fitness_ = 1.0;

    //! Extra data
    size_t id_ = 0;
    size_t mother_ = 0;
    size_t ancestor_ = 0;
    double time_of_birth_ = 0.0;
    double time_of_death_ = 0.0;
};

#endif /* GLAND_H_ */
