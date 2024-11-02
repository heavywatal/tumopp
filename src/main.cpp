/*! @file main.cpp
    @brief Defines main() and file output
*/
#include "simulation.hpp"
#include <iostream>

//! Instantiate and run Simulation
int main(int argc, char* argv[]) {
    std::vector<std::string> arguments(argv + 1, argv + argc);
    try {
        tumopp::Simulation simulation(arguments);
        simulation.run();
        simulation.write();
    } catch (const tumopp::exit_success& e) {
        return 0;
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }
    return 0;
}
