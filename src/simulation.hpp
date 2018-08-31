/*! @file simulation.hpp
    @brief Interface of Simulation class
*/
#pragma once
#ifndef SIMULATION_HPP_
#define SIMULATION_HPP_

#include <iosfwd>
#include <vector>
#include <string>
#include <memory>

namespace boost {namespace program_options {
  class options_description;
  class variables_map;
}}

namespace tumopp {

class Tissue;
struct EventRates;
struct CellParams;

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

    //! Write files unless `--outdir` is empty
    void write() const;

    //! @name Output for Rcpp
    //@{
    std::string history() const;
    std::string snapshots() const;
    std::string drivers() const;
    std::string config_string() const {return config_string_;}
    //@}
    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
  private:
    boost::program_options::options_description options_desc();
    boost::program_options::options_description cell_options();

    //! Print help message and exit
    [[noreturn]] void help_and_exit();
    //! Print version and exit
    [[noreturn]] void version_and_exit();

    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
    // Data member

    //! Tissue instance
    std::unique_ptr<Tissue> tissue_;
    //! EventRates instance
    std::unique_ptr<EventRates> init_event_rates_;
    //! CellParams instance
    std::unique_ptr<CellParams> cell_params_;

    //! optional variables
    std::unique_ptr<boost::program_options::variables_map> vars_;
    //! Command line arguments
    std::string command_args_;
    //! Writen to "program_options.conf"
    std::string config_string_;
};

} // namespace tumopp

#endif /* SIMULATION_HPP_ */
