/*! @file simulation.hpp
    @brief Interface of Simulation class
*/
#pragma once
#ifndef TUMOPP_SIMULATION_HPP_
#define TUMOPP_SIMULATION_HPP_

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

    //! @name Getter for `main()`
    //@{
    Tissue& tissue() const noexcept {return *tissue_;}
    const std::string& config() const noexcept {return config_;}
    std::string outdir() const;
    //@}

    //! @name Output for Rcpp
    //@{
    std::string history() const;
    std::string snapshots() const;
    std::string drivers() const;
    std::string benchmark() const;
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

#endif // TUMOPP_SIMULATION_HPP_
