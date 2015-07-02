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
    static size_t MAX_SIZE_;
    static size_t DIMENSIONS_;

    Tissue() {
        tumor_.emplace(std::vector<int>(DIMENSIONS_), Gland());
    }

    std::string grow();

    std::string tsv() const;

    static void unit_test();
    static boost::program_options::options_description& opt_description();

  private:
    void push(const Gland& daughter, std::vector<int> parent_pos, const std::vector<int>& direction);

    std::map<std::vector<int>, Gland> tumor_;
};

#endif /* TISSUE_H_ */
