// -*- mode: c++; coding: utf-8 -*-
/*! @file tissue.h
    @brief Interface of Tissue class
*/
#pragma once
#ifndef TISSUE_H_
#define TISSUE_H_
#include <cmath>
#include <iostream>
#include <vector>
#include <map>
#include <string>

#include "gland.h"

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

namespace boost {
    namespace program_options {
        class options_description;
    }
}

class Tissue {
  public:
    //! Constructor
    Tissue(const std::vector<int>& origin): coords_{1, origin} {
        tumor_.emplace(origin, Gland());
        init_regularly();
    }
    //! Constructor
    Tissue(const size_t dimensions=2): Tissue{std::vector<int>(dimensions)} {};

    //! Initiate tumor with 2^D glands regularly
    void init_regularly();
    //! Initiate tumor with 2^D glands randomly
    void init_randomly();
    //! Mark first cells with mutation
    void stain();

    //! Big bang model
    void grow_random(const size_t max_size);
    //! Another extreme
    void grow_even(const size_t max_size);

    //! Return tumor state as TSV string
    std::string snapshot(const std::string& sep="\t") const;
    //! Return mutation_coords_ as TSV string
    std::string mutation_history(const std::string& sep="\t") const;

    friend std::ostream& operator<< (std::ostream&, const Tissue&);

    //! Unit test
    static void unit_test();
    static boost::program_options::options_description& opt_description();

  private:
    //! Emplace daughter gland at the specified coord
    void emplace(const std::vector<int>& current_coord, Gland&& daughter);
    //! Emplace daughter gland and push other glands outward
    void push(Gland&& daughter, std::vector<int>* current_coord, const std::vector<int>& direction);
    //! Emplace daughter gland at a neighboring site and push if needed
    void fill_neighbor(Gland&& daughter, const std::vector<int>& current_coord);

    //! key: coords, value: gland
    std::map<std::vector<int>, Gland> tumor_;
    //! The coordinates of the existing glands
    std::vector<std::vector<int>> coords_;
    //! The coordinates of the past mutations
    std::vector<std::vector<int>> mutation_coords_;
    //! Timing of mutations (tumor size as proxy)
    std::vector<size_t> mutation_stages_;
};

#endif /* TISSUE_H_ */
