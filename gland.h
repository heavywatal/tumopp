// -*- mode: c++; coding: utf-8 -*-
/*! @file gland.h
    @brief Interface of Gland class
*/
#pragma once
#ifndef GLAND_H_
#define GLAND_H_
#include <iostream>
#include <vector>
#include <bitset>
#include <string>

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

namespace boost {
    namespace program_options {
        class options_description;
    }
}

class Gland {
  public:
    static double CELLS_PER_GLAND_;
    static double MUTATION_RATE_;
    static double MUTATION_SIGMA_;
    static double APOPTOSIS_RATE_;

    bool mutate();
    bool apoptosis() const;

    double fitness() const {return fitness_;}
    const std::vector<size_t>& sites() const {return sites_;}

    friend std::ostream& operator<< (std::ostream& ost, const Gland& gland) {
        return ost << gland.fitness();
    }

    static void unit_test();
    static boost::program_options::options_description& opt_description();

  private:
    static size_t MUTATED_SITES_;

    double fitness_ = 1.0;
    std::vector<size_t> sites_;
};

#endif /* GLAND_H_ */
