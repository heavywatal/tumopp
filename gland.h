// -*- mode: c++; coding: utf-8 -*-
/*! @file gland.h
    @brief Interface of Gland class
*/
#pragma once
#ifndef GLAND_H_
#define GLAND_H_

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
    //! Constructor
    Gland(const std::vector<int>& v): coord_(v), id_(++ID_TAIL_), ancestor_(id_) {}
    //! Copy constructor
    Gland(const Gland& other):
        coord_(other.coord_), sites_(other.sites_), age_(other.age_),
        id_(++ID_TAIL_), mother_(other.id_), ancestor_(other.ancestor_) {}
    //! Copy assignment operator
    Gland& operator=(const Gland&) = delete;
    //! Move constructor
    Gland(Gland&& other) = default;
    //! Move assignment operator
    Gland& operator=(Gland&&) = default;

    void set_coord(const std::vector<int>& v) {coord_ = v;}

    //! Calculate fitness
    double fitness() const {
        double result = 1.0;
        for (const auto i: sites_) {
            result += MUTATION_EFFECTS_[i];
        }
        return result;
        //return std::min(std::max(1.0, result), 5.0);
    }

    //! Mutate and record
    void mutate();

    //! Update age_
    void stamp(const size_t x) {age_ = x;}

    //! Bernoulli trial of mutation
    static bool bernoulli_mutation();
    //! Bernoulli trial of birth
    bool bernoulli_birth() const;
    //! Bernoulli trial of death
    bool bernoulli_death() const;

    //! Getter
    const std::vector<int>& coord() const {return coord_;}
    //! Getter
    const std::vector<size_t>& sites() const {return sites_;}
    //! Getter
    size_t age() const {return age_;}
    size_t id() const {return id_;}
    size_t mother() const {return mother_;}
    size_t ancestor() const {return ancestor_;}
    //! Getter
    static const std::vector<double>& MUTATION_EFFECTS() {return MUTATION_EFFECTS_;}

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

    //! Position in a tumor
    std::vector<int> coord_;
    //! Mutated sites (infinite-site model)
    std::vector<size_t> sites_;
    //! The age of the last division
    size_t age_ = 0;
    size_t id_ = 0;
    size_t mother_ = 0;
    size_t ancestor_ = 0;
};

#endif /* GLAND_H_ */
