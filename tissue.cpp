// -*- mode: c++; coding: utf-8 -*-
/*! @file tissue.cpp
    @brief Implementation of Gland class
*/
#include "tissue.h"

#include <cmath>
#include <iostream>
#include <sstream>

#include <boost/program_options.hpp>

#include "cxxwtils/prandom.hpp"
#include "cxxwtils/iostr.hpp"

size_t Tissue::MAX_SIZE_ = 20;

boost::program_options::options_description& Tissue::opt_description() {
    namespace po = boost::program_options;
    static po::options_description desc{"Tissue"};
    desc.add_options()
        ("max,l", po::value<size_t>(&MAX_SIZE_)->default_value(MAX_SIZE_))
    ;
    return desc;
}

void Tissue::grow() {
    glands_.resize(MAX_SIZE_);
    glands_[MAX_SIZE_ / 2 - 1] = Gland(true);
    std::cout << str_malignancy() << std::endl;
    for (size_t t=0; t<20; ++t) {
        for (size_t i=0; i<glands_.size(); ++i) {
            auto& parent = glands_[i];
            if (!parent.malignant()) continue;
            Gland daughter = parent;
            daughter.mutate();
            if (parent.apoptosis()) {
                parent = daughter;
            } else {
                const size_t direction = static_cast<size_t>(wtl::prandom().bernoulli(0.5));
                auto pos = i + direction;
                if (pos < MAX_SIZE_ / 2) {
                    glands_.pop_front();
                } else {
                    glands_.pop_back();
                    ++i;
                }
                glands_.insert(glands_.begin() + pos, daughter);
            }
        }
        //std::cout << str_fitness() << std::endl;
        std::cout << str_malignancy() << std::endl;
    }
}

std::string Tissue::str_fitness() const {
    return wtl::oss_join(glands_, ",", [](const Gland& x) {return x.fitness();});
}

std::string Tissue::str_malignancy() const {
    return wtl::oss_join(glands_, "", [](const Gland& x) {return x.malignant();});
}

void Tissue::unit_test() {
    std::cerr << __PRETTY_FUNCTION__ << std::endl;
    std::cerr.precision(15);
}
