// -*- mode: c++; coding: utf-8 -*-
/*! @file simulation.h
    @brief Interface of Simulation class
*/
#pragma once
#ifndef SIMULATION_H_
#define SIMULATION_H_

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

/*! @brief Represents single run
*/
class Simulation {
  public:
    //! Parse command arguments
    Simulation(int argc, char* argv[]);

    //! Top level function that should be called from main()
    void run();

    boost::program_options::options_description& opt_description();

    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
  private:

    //! Max tumor size to stop simulation
    size_t MAX_SIZE = 16384;

    //! Coordinate system
    std::string COORDINATE = "moore";

    //! Print extra information
    bool VERBOSE = false;

    //! Seed for random number generator
    unsigned int SEED = std::random_device{}();

    //! Control execution mode
    int MODE = 0;

    //! Group name of this run such as altered parameter
    std::string LABEL;

    //! Home directory (c_str)
    const char* HOME_CHAR = std::getenv("HOME");

    //! Home directory
    const fs::path HOME_DIR = HOME_CHAR;

    //! Temporary directory
    const fs::path TMP_DIR = HOME_DIR / "tmp";

    //! Working directory to write out temporal results
    fs::path WORK_DIR;

    //! Target directory to which the contents in WORK_DIR are moved
    fs::path OUT_DIR = TMP_DIR / wtl::strftime("crc%Y%m%d");

    //! Project directory
    fs::path PROJECT_DIR;
};

#endif /* SIMULATION_H_ */
