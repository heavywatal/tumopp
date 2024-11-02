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
#include <stdexcept>

namespace tumopp {

class Tissue;
struct EventRates;
struct CellParams;

class exit_success: public std::logic_error {
  public:
    exit_success() noexcept: std::logic_error("") {}
};

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

    //! Write results to files
    void write() const;

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

} // namespace tumopp

#endif // TUMOPP_SIMULATION_HPP_
