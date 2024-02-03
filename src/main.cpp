/*! @file main.cpp
    @brief Defines main() and file output
*/
#include "simulation.hpp"
#include "tissue.hpp"

#include <wtl/filesystem.hpp>
#ifdef ZLIB_FOUND
  #include <wtl/zlib.hpp>
#endif

#include <iostream>
#include <fstream>

//! Write config and simulation result to files
void write(tumopp::Simulation& simulation) {
  #ifdef ZLIB_FOUND
    using ofstream = wtl::zlib::ofstream;
    const std::string ext = ".tsv.gz";
  #else
    using ofstream = std::ofstream;
    const std::string ext = ".tsv";
  #endif
    namespace fs = wtl::filesystem;
    const auto outdir = simulation.outdir();
    if (outdir.empty()) return;
    fs::create_directory(outdir);
    fs::current_path(outdir);
    std::ofstream{"config.json"} << simulation.config();
    const auto& tissue = simulation.tissue();
    {
        ofstream ofs{"population" + ext};
        ofs.exceptions(std::ios_base::failbit | std::ios_base::badbit);
        tissue.write_history(ofs);
    }
    if (tissue.has_snapshots()) {
        ofstream ofs{"snapshots" + ext};
        ofs.exceptions(std::ios_base::failbit | std::ios_base::badbit);
        tissue.write_snapshots(ofs);
    }
    if (tissue.has_drivers()) {
        ofstream ofs{"drivers" + ext};
        ofs.exceptions(std::ios_base::failbit | std::ios_base::badbit);
        tissue.write_drivers(ofs);
    }
    if (tissue.has_benchmark()) {
        ofstream ofs{"benchmark" + ext};
        ofs.exceptions(std::ios_base::failbit | std::ios_base::badbit);
        tissue.write_benchmark(ofs);
    }
}

//! Instantiate and run Simulation
int main(int argc, char* argv[]) {
    std::vector<std::string> arguments(argv + 1, argv + argc);
    try {
        tumopp::Simulation simulation(arguments);
        simulation.run();
        write(simulation);
    } catch (const std::runtime_error& e) {
        std::cerr << e.what() << std::endl;
    }
    return 0;
}
