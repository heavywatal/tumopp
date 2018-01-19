/*! @file simulation.hpp
    @brief Interface of Simulation class
*/
#pragma once
#ifndef SIMULATION_HPP_
#define SIMULATION_HPP_

#include "tissue.hpp"

#include <wtl/exception.hpp>
#include <boost/program_options.hpp>

#include <iostream>
#include <sstream>
#include <vector>
#include <random>

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

    //! Write ms-like sequence and files (if #WRITE_TO_FILES)
    void write() const;

    //! getter of #tissue_ for tumorr
    const Tissue& tissue() const {return tissue_;}

    //! getter of #CONFIG_STRING for tumorr
    const std::string& conf() const {return CONFIG_STRING;}

    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
  private:
    //! options description for optional arguments
    boost::program_options::options_description options_desc();
    //! options description for positional arguments
    boost::program_options::options_description positional_desc();
    //! Print help message and exit
    [[noreturn]] void help_and_exit();

    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
    // Data member

    //! Tissue instance
    Tissue tissue_;

    //! Number of samples
    size_t NSAM = 20;
    //! Number of repeats
    size_t HOWMANY = 1;
    //! Max tumor size to stop simulation
    size_t MAX_SIZE = 16384;
    //! Duration of turnover phase after log growth
    double PLATEAU_TIME = 0.0;
    //! -w
    bool WRITE_TO_FILES = false;
    //! Seed for random number generator
    unsigned int SEED = std::random_device{}();

    //! Output directory is created if #WRITE_TO_FILES
    std::string OUT_DIR;
    //! Command line arguments
    std::string COMMAND_ARGS;
    //! Writen to "program_options.conf"
    std::string CONFIG_STRING;
};

} // namespace tumopp

#endif /* SIMULATION_HPP_ */
