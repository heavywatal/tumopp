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

#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

namespace boost {
    namespace program_options {
        class options_description;
    }
}

class Tissue;

/*! @brief Represents single run
*/
class Simulation {
  public:
    //! Parse command arguments
    Simulation(int argc, char* argv[]);

    //! Top level function that should be called once from main()
    void run() const;

    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
  private:
    boost::program_options::options_description& opt_description();

    //! Max tumor size to stop simulation
    size_t MAX_SIZE = 16384;

    size_t NSAM = 20;

    size_t HOWMANY = 1;

    //! Print extra information
    bool VERBOSE = false;

    //! Seed for random number generator
    unsigned int SEED = std::random_device{}();

    //! Target directory to which the contents in WORK_DIR are moved
    fs::path OUT_DIR;

    std::string COMMAND_ARGS;
    std::string CONFIG_STRING;
};

#endif /* SIMULATION_HPP_ */