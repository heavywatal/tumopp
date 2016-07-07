// -*- mode: c++; coding: utf-8 -*-
/*! @file simulation.cpp
    @brief Inplementation of Simulation class
*/
#include "simulation.hpp"

#include <cstdlib>
#include <algorithm>

#include <boost/filesystem.hpp>

#include <cxxwtils/iostr.hpp>
#include <cxxwtils/getopt.hpp>
#include <cxxwtils/prandom.hpp>
#include <cxxwtils/os.hpp>
#include <cxxwtils/gz.hpp>

#include "cell.hpp"

namespace fs = boost::filesystem;

namespace tumopp {

boost::program_options::options_description Simulation::opt_description() {HERE;
    namespace po = boost::program_options;
    po::options_description description("Simulation");
    description.add_options()
        ("help,h", po::value<bool>()->default_value(false)->implicit_value(true), "produce help")
        ("verbose,v", po::value<bool>(&VERBOSE)
            ->default_value(VERBOSE)->implicit_value(true), "verbose output")
        ("test", po::value<int>()->default_value(0)->implicit_value(1))
        ("out_dir,o", po::value<std::string>(&OUT_DIR)->default_value(OUT_DIR))
        ("seed", po::value<unsigned int>(&SEED)->default_value(SEED))
        ("nsam", po::value<size_t>(&NSAM)->default_value(NSAM))
        ("howmany", po::value<size_t>(&HOWMANY)->default_value(HOWMANY))
    ;
    return description;
}

//! Unit test for each class
inline void test(const int flg) {HERE;
    switch (flg) {
      case 0:
        break;
      case 1:
        Cell::unit_test();
        Tissue::unit_test();
        std::exit(EXIT_SUCCESS);
      default:
        std::abort();
    }
}

Simulation::Simulation(const std::vector<std::string>& arguments) {HERE;
    COMMAND_ARGS = wtl::str_join(arguments, " ") ;
    OUT_DIR = wtl::strftime("tumopp_%Y%m%d_%H%M_") + std::to_string(::getpid());

    namespace po = boost::program_options;
    po::options_description description;
    description.add(opt_description());
    description.add(Cell::opt_description());
    description.add(Tissue::opt_description());

    po::positional_options_description positional;
    positional.add("nsam", 1);
    positional.add("howmany", 1);

    po::variables_map vm;
    po::store(po::command_line_parser(arguments).
              options(description).
              positional(positional).run(), vm);
    po::notify(vm);
    Cell::init_distributions();

    if (vm["help"].as<bool>()) {
        description.print(std::cout);
        std::exit(EXIT_SUCCESS);
    }
    wtl::sfmt().seed(SEED);
    CONFIG_STRING = wtl::flags_into_string(vm);
    if (VERBOSE) {
        std::cerr << wtl::iso8601datetime() << std::endl;
        std::cerr << CONFIG_STRING << std::endl;
    }
    test(vm["test"].as<int>());
    if (NSAM > vm["max"].as<size_t>()) {
        std::cout.precision(1);
        std::cout
            << "\nERROR: NSAM=" << NSAM
            << " is larger than tumor size "
            << std::fixed << vm["max"].as<size_t>() << std::endl;
        std::abort();
    }
//    OUT_DIR = fs::path(vm["out_dir"].as<std::string>());
    OUT_DIR = fs::system_complete(OUT_DIR).string();
}

void Simulation::run() {HERE;
    std::ios::sync_with_stdio(false);
    std::cin.tie(0);
    std::cout.precision(15);
    std::cerr.precision(6);
    while (!tissue_.grow()) {;}
}

void Simulation::write() const {HERE;
    std::cout << "tumopp " << COMMAND_ARGS << "\n" << SEED << "\n";

    if (Tissue::DIMENSIONS() == 3) {
        for (size_t i=0; i<HOWMANY; ++i) {
            tissue_.write_segsites(std::cout, tissue_.sample_section(NSAM));
        }
    } else {
        for (size_t i=0; i<HOWMANY; ++i) {
            tissue_.write_segsites(std::cout, tissue_.sample_random(NSAM));
        }
    }

    if (VERBOSE) {
        derr("mkdir && cd to " << OUT_DIR << std::endl);
        fs::create_directory(OUT_DIR);
        wtl::cd(OUT_DIR);
        wtl::Fout{"program_options.conf"} << CONFIG_STRING;
        wtl::gzip{wtl::Fout{"population.tsv.gz"}}
            << tissue_.specimens();
        wtl::gzip{wtl::Fout{"snapshots.tsv.gz"}}
            << tissue_.snapshots();
        wtl::gzip{wtl::Fout{"drivers.tsv.gz"}}
            << tissue_.drivers();
        wtl::gzip{wtl::Fout{"distance.tsv.gz"}}
            << tissue_.pairwise_distance(std::min(200UL, tissue_.size()));
        std::cerr << wtl::iso8601datetime() << std::endl;
    }
}

} // namespace tumopp
