/*! @file simulation.hpp
    @brief Interface of Simulation class
*/
#pragma once
#ifndef SIMULATION_HPP_
#define SIMULATION_HPP_

#include "tissue.hpp"

#include <vector>
#include <string>
#include <memory>

namespace boost {namespace program_options {
  class options_description;
  class variables_map;
}}

namespace tumopp {

class Tissue;

/*! @brief Represents single run
*/
class Simulation {
  public:
    //! Parse command arguments
    Simulation(const std::vector<std::string>& args);

    //! Non-default destructor for forward declaration
    ~Simulation();

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

    //! optional variables
    std::unique_ptr<boost::program_options::variables_map> vars_;
    //! Command line arguments
    std::string command_args_;
    //! Writen to "program_options.conf"
    std::string config_string_;
};

} // namespace tumopp

#endif /* SIMULATION_HPP_ */
