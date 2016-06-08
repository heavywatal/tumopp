// -*- mode: c++; coding: utf-8 -*-
/*! @file main.cpp
    @brief Only defines tiny main()
*/
#include "simulation.hpp"

//! Just instantiate and run Simulation
int main(int argc, char* argv[]) {
    Simulation simulation(argc, argv);
    simulation.run();
    return 0;
}

