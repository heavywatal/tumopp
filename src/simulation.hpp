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

    //! Make samples and write ms-like output
    void ms(std::ostream&) const;

    //! Shortcut of Tissue.specimens() for R
    std::string specimens() const;
    //! Shortcut of Tissue.snapshots() for R
    std::string snapshots() const;
    //! Shortcut of Tissue.drivers() for R
    std::string drivers() const;
    //! Shortcut of Tissue.pairwise_distance() for R
    std::string pairwise_distance(size_t npair) const;
    //! Shortcut of ms() as std::string for R
    std::string ms() const;
    //! getter of #config_string_
    std::string config_string() const {return config_string_;}

    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
  private:
    boost::program_options::options_description options_desc();
    boost::program_options::options_description positional_desc();

    //! Print help message and exit
    [[noreturn]] void help_and_exit();
    //! Print version and exit
    [[noreturn]] void version_and_exit();

    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
    // Data member

    //! Tissue instance
    std::unique_ptr<Tissue> tissue_;

    //! optional variables
    std::unique_ptr<boost::program_options::variables_map> vars_;
    //! Command line arguments
    std::string command_args_;
    //! Writen to "program_options.conf"
    std::string config_string_;
};

} // namespace tumopp

#endif /* SIMULATION_HPP_ */
