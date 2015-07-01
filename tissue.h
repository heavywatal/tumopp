// -*- mode: c++; coding: utf-8 -*-
/*! @file tissue.h
    @brief Interface of Gland class
*/
#pragma once
#ifndef TISSUE_H_
#define TISSUE_H_
#include <cmath>
#include <iostream>
#include <vector>
#include <deque>
#include <bitset>
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
    static size_t MAX_SIZE_;

    void grow();

    static void unit_test();
    static boost::program_options::options_description& opt_description();

    std::string str_fitness() const;
    std::string str_malignancy() const;

  private:
    std::deque<Gland> glands_;
};

#endif /* TISSUE_H_ */
