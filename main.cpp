/*! @file main.cpp
    @brief Only defines tiny main()
*/
#include "src/simulation.hpp"

//! Just instantiate and run Simulation
int main(int argc, char* argv[]) {
    std::vector<std::string> arguments(argv + 1, argv + argc);
    tumopp::Simulation simulation(arguments);
    simulation.run();
    simulation.write();
    return 0;
}
