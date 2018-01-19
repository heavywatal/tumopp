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

    //! Write ms-like sequence and files
    void write() const;

    //! getter of #tissue_ for tumorr
    const Tissue& tissue() const {return tissue_;}

    //! getter of #config_string_ for tumorr
    const std::string& conf() const {return config_string_;}

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
    size_t nsam_ = 20;
    //! Number of repeats
    size_t howmany_ = 1;
    //! Max tumor size to stop simulation
    size_t max_size_ = 16384;
    //! Duration of turnover phase after log growth
    double plateau_time_ = 0.0;
    //! Seed for random number generator
    unsigned int seed_ = std::random_device{}();

    //! Output directory
    std::string out_dir_;
    //! Command line arguments
    std::string command_args_;
    //! Writen to "program_options.conf"
    std::string config_string_;
};

} // namespace tumopp

#endif /* SIMULATION_HPP_ */
