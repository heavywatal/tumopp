/*! @file main.cpp
    @brief Only defines tiny main()
*/
#include "src/simulation.hpp"
#include <iostream>

//! Just instantiate and run Simulation
int main(int argc, char* argv[]) {
    std::vector<std::string> arguments(argv + 1, argv + argc);
    try {
        tumopp::Simulation simulation(arguments);
        simulation.run();
        simulation.write();
    } catch (const std::runtime_error& e) {
        std::cerr << e.what() << std::endl;
    }
    return 0;
}
