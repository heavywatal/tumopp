// -*- mode: c++; coding: utf-8 -*-
/*! @file cell.h
    @brief Interface of Cell class
*/
#pragma once
#ifndef CELL_HPP_
#define CELL_HPP_

#include <cmath>
#include <vector>
#include <string>

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

namespace boost {
    namespace program_options {
        class options_description;
    }
}

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
    Cell(const std::vector<int>& v): coord_(v), id_(++ID_TAIL_) {}
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
    void mutate();

    //! Setter
    void set_coord(const std::vector<int>& v) {coord_ = v;}
    void set_time_of_birth(const double t) {
        time_of_birth_ = t;
        id_ = ++ID_TAIL_;
        if (type_ == CellType::nonstem) {--proliferation_capacity_;}
    }
    void set_time_of_death(const double t) {time_of_death_ = t;}
    void daughterize(const double t) {
        time_of_death_ = 0.0;
        ancestors_.push_back(id_);
        set_time_of_birth(t);
    }

    //! Getter
    bool is_dividing() const {return next_event_ == Event::birth;}
    bool is_dying() const {return next_event_ == Event::death;}
    bool is_migrating() const {return next_event_ == Event::migration;}
    double delta_time(const double positional_value);

    const std::vector<int>& coord() const {return coord_;}
    const std::vector<size_t>& sites() const {return sites_;}
    size_t id() const {return id_;}
    double time_of_birth() const {return time_of_birth_;}
    double time_of_death() const {return time_of_death_;}

    std::vector<int> is_descendant_of(const std::vector<size_t>&);

    //! convert site positions to 01 vector
    std::vector<size_t> haplotype(std::vector<size_t> segsites) const;

    static std::vector<size_t> GENERATE_NEUTRAL_MUTATIONS();

    static std::string header(const size_t dimensions, const char* sep);
    std::ostream& write(std::ostream& ost, const char* sep) const;
    std::string str(const char* sep) const;
    friend std::ostream& operator<< (std::ostream&, const Cell&);

    //! Unit test
    static void unit_test();
    static boost::program_options::options_description& opt_description();

  private:
    //! finite per capita rates
    double growth_rate() const {return 1.0 + birth_rate_ - death_rate_;}
    //! instantaneous rate for time increment
    double instantaneous_event_rate() const {
        const double lambda = growth_rate();
        return (birth_rate_ + death_rate_) * std::log(lambda) / (lambda - 1.0);
    }

    //! per cell division
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

    static double BIRTH_RATE_;

    static double DEATH_RATE_;

    static double MIGRATION_RATE_;
    //! k
    static double GAMMA_SHAPE_;

    static double PROB_SYMMETRIC_DIVISION_;
    static size_t MAX_PROLIFERATION_CAPACITY_;

    static size_t ID_TAIL_;

    //! Position in a tumor
    std::vector<int> coord_;
    //! Mutated sites (infinite-site model)
    std::vector<size_t> sites_;

    double birth_rate_ = BIRTH_RATE_;
    double death_rate_ = DEATH_RATE_;
    double migra_rate_ = MIGRATION_RATE_;

    CellType type_ = CellType::stem;
    size_t proliferation_capacity_ = MAX_PROLIFERATION_CAPACITY_;

    //! Extra data
    size_t id_ = 0;
    std::vector<size_t> ancestors_;
    double time_of_birth_ = 0.0;
    double time_of_death_ = 0.0;
    Event next_event_ = Event::birth;
};

#endif /* CELL_HPP_ */
