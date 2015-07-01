// -*- mode: c++; coding: utf-8 -*-
/*! @file main.cpp
    @brief Only defines tiny main()
*/
#include <iostream>

#include "simulation.h"

//! Just instantiate and run Simulation
int main(int argc, char* argv[]) {
    std::ios::sync_with_stdio(false);
    std::cin.tie(0);
    std::cout.precision(16);
    std::cerr.precision(6);

    Simulation simulation(argc, argv);
    simulation.run();
    return 0;
}

