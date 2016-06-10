// -*- mode: c++; coding: utf-8 -*-
/*! @file main.cpp
    @brief Only defines tiny main()
*/
#include <tumopp/simulation.hpp>

//! Just instantiate and run Simulation
int main(int argc, char* argv[]) {
    tumopp::Simulation simulation(argc, argv);
    simulation.run();
    return 0;
}

