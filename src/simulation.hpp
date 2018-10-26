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
    std::string benchmark() const;
    std::string config() const {return config_;}
    //@}
    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
  private:
    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
    // Data member

    //! Tissue instance
    std::unique_ptr<Tissue> tissue_;
    //! EventRates instance
    std::unique_ptr<EventRates> init_event_rates_;
    //! CellParams instance
    std::unique_ptr<CellParams> cell_params_;
    //! Parameters
    std::string config_;
};

//! @name Workaround for R/Rcpp
//@{
std::streambuf* std_cout_rdbuf(std::streambuf*);
std::streambuf* std_cerr_rdbuf(std::streambuf*);
//@}

} // namespace tumopp

#endif /* SIMULATION_HPP_ */
