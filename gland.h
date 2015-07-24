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
    //! 80 billion cells / 8 million glands = 10000
    static double CELLS_PER_GLAND_;

    //! per gland per division {1e-8 to 1e-4}
    static double MUTATION_RATE_;

    //! {0, 0.2, 0.6}
    static double MUTATION_SIGMA_;

    //! 0.2 fixed in the article
    static double APOPTOSIS_RATE_;

    //! Default constructor
    Gland() = default;
    //! Copy constructor
    Gland(const Gland& other) = default;
    //! Copy assignment operator
    Gland& operator=(const Gland&) = delete;
    //! Move constructor
    Gland(Gland&& other) = default;
    //! Move assignment operator
    Gland& operator=(Gland&&) = default;

    //! Calculate fitness
    double fitness() const {
        double result = 1.0;
        for (const auto i: sites_) {
            result += MUTATION_EFFECTS_[i];
        }
        return std::min(std::max(1.0, result), 5.0);
    }

    //! Mutate and record
    void mutate();

    //! Update age_
    void stamp(const size_t x) {age_ = x;}

    //! Bernoulli trial of mutation
    static bool bernoulli_mutation();
    //! Bernoulli trial of apoptosis
    bool bernoulli_apoptosis() const;

    //! Getter
    const std::vector<size_t>& sites() const {return sites_;}
    //! Getter
    size_t age() const {return age_;}
    //! Getter
    static const std::vector<double>& MUTATION_EFFECTS() {return MUTATION_EFFECTS_;}

    friend std::ostream& operator<< (std::ostream&, const Gland&);

    //! Unit test
    static void unit_test();
    static boost::program_options::options_description& opt_description();

  private:
    //! The history of mutation effects
    static std::vector<double> MUTATION_EFFECTS_;

    //! Mutated sites (infinite-site model)
    std::vector<size_t> sites_;
    //! The age of the last division
    size_t age_ = 0;
};

#endif /* GLAND_H_ */
