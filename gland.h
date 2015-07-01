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
    static double MUTATION_RATE_;
    static double SIGMA_;

    Gland(bool malignant=false): malignant_{malignant} {}

    void mutate();
    bool apoptosis() const;
    int fission() const;
    double fitness() const {return fitness_;}
    unsigned int malignant() const {return static_cast<unsigned int>(malignant_);}

    static void unit_test();
    static boost::program_options::options_description& opt_description();

  private:
    static size_t MUTATED_SITES_;

    bool malignant_ = false;
    double fitness_ = 1.0;
    std::vector<size_t> sites_;
};

#endif /* GLAND_H_ */
