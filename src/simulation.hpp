// -*- mode: c++; coding: utf-8 -*-
/*! @file simulation.h
    @brief Interface of Simulation class
*/
#pragma once
#ifndef SIMULATION_HPP_
#define SIMULATION_HPP_

#include <iostream>
#include <sstream>
#include <vector>
#include <random>

#include <cxxwtils/iostr.hpp>

#include "tissue.hpp"

namespace boost {
    namespace program_options {
        class options_description;
    }
}

namespace tumopp {

class Tissue;

/*! @brief Represents single run
*/
class Simulation {
  public:
    //! Parse command arguments
    Simulation(const std::vector<std::string>& args);

    //! Top level function that should be called once from main()
    void run();

    void write() const;

    const Tissue& tissue() const {return tissue_;}

    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
  private:
    boost::program_options::options_description& opt_description();

    Tissue tissue_;

    size_t NSAM = 20;

    size_t HOWMANY = 1;

    //! Print extra information
    bool VERBOSE = false;

    //! Seed for random number generator
    unsigned int SEED = std::random_device{}();

    //! Target directory to which the contents in WORK_DIR are moved
    std::string OUT_DIR;

    std::string COMMAND_ARGS;
    std::string CONFIG_STRING;
};

} // namespace tumopp

#endif /* SIMULATION_HPP_ */
