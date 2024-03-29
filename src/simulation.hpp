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
    //! Get #tissue_
    Tissue& tissue() const noexcept {return *tissue_;}
    //! Get #config_
    const std::string& config() const noexcept {return config_;}
    //! Get VM["outdir"]
    std::string outdir() const;
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

} // namespace tumopp

#endif // TUMOPP_SIMULATION_HPP_
