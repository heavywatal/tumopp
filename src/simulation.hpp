// -*- mode: c++; coding: utf-8 -*-
/*! @file simulation.hpp
    @brief Interface of Simulation class
*/
#pragma once
#ifndef SIMULATION_HPP_
#define SIMULATION_HPP_

#include <iostream>
#include <sstream>
#include <vector>
#include <random>

#include <wtl/exception.hpp>

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

    std::string conf() const {return CONFIG_STRING;}

    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
  private:
    boost::program_options::options_description options_desc();
    boost::program_options::options_description positional_desc();
    void help_and_exit();

    Tissue tissue_;

    size_t NSAM = 20;

    size_t HOWMANY = 1;

    //! Max tumor size to stop simulation
    size_t MAX_SIZE = 16384;

    bool WRITE_TO_FILES = false;

    //! Seed for random number generator
    unsigned int SEED = std::random_device{}();

    //! Target directory to which the contents in WORK_DIR are moved
    std::string OUT_DIR;

    std::string COMMAND_ARGS;
    std::string CONFIG_STRING;
};

} // namespace tumopp

#endif /* SIMULATION_HPP_ */
