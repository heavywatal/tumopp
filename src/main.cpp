/*! @file main.cpp
    @brief Only defines tiny main()
*/
#include "simulation.hpp"
#include "tissue.hpp"

#include <wtl/debug.hpp>
#include <wtl/filesystem.hpp>
#include <wtl/zlib.hpp>

#include <iostream>

//! variables mapper to handle command-line arguments
namespace fs = wtl::filesystem;

void write(tumopp::Simulation& simulation) {
    const auto& tissue = simulation.tissue();
    const auto outdir = simulation.outdir();
    if (outdir.empty()) return;
    DCERR("mkdir && cd to " << outdir << std::endl);
    fs::create_directory(outdir);
    fs::current_path(outdir);
    std::cerr << "Output: " << outdir << "\n";
    std::ofstream{"config.json"} << simulation.config();
    {
        wtl::zlib::ofstream ofs{"population.tsv.gz"};
        tissue.write_history(ofs);
    }
    if (tissue.has_snapshots()) {
        wtl::zlib::ofstream ofs{"snapshots.tsv.gz"};
        tissue.write_snapshots(ofs);
    }
    if (tissue.has_drivers()) {
        wtl::zlib::ofstream ofs{"drivers.tsv.gz"};
        tissue.write_drivers(ofs);
    }
    if (tissue.has_benchmark()) {
        wtl::zlib::ofstream ofs{"benchmark.tsv.gz"};
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
