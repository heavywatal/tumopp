/*! @file main.cpp
    @brief Only defines tiny main()
*/
#include "simulation.hpp"
#include "tissue.hpp"

#include <wtl/filesystem.hpp>
#ifdef ZLIB_FOUND
  #include <wtl/zlib.hpp>
#endif

#include <iostream>
#include <fstream>

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
    std::cerr << "Output: " << fs::current_path() << "\n";
    std::ofstream{"config.json"} << simulation.config();
    const auto& tissue = simulation.tissue();
    {
        ofstream ofs{"population" + ext};
        tissue.write_history(ofs);
    }
    if (tissue.has_snapshots()) {
        ofstream ofs{"snapshots" + ext};
        tissue.write_snapshots(ofs);
    }
    if (tissue.has_drivers()) {
        ofstream ofs{"drivers" + ext};
        tissue.write_drivers(ofs);
    }
    if (tissue.has_benchmark()) {
        ofstream ofs{"benchmark" + ext};
        tissue.write_benchmark(ofs);
    }
}

//! Just instantiate and run Simulation
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
