// -*- mode: c++; coding: utf-8 -*-
/*! @file simulation.cpp
    @brief Inplementation of Simulation class
*/
#include "simulation.h"

#include <cxxwtils/iostr.hpp>
#include <cxxwtils/getopt.hpp>
#include <cxxwtils/prandom.hpp>
#include <cxxwtils/os.hpp>
#include <cxxwtils/gz.hpp>

#include "tissue.h"
#include "gland.h"

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
// functions

//! Program options
/*! @return Program options description

    Command line option | Symbol | Variable
    --------------------| ------ | -------------------------
    `-N,--max`          | -      | Simulation::MAX_SIZE
*/
boost::program_options::options_description& Simulation::opt_description() {HERE;
    namespace po = boost::program_options;
    static po::options_description description("Simulation");
    description.add_options()
        ("max,N", po::value<size_t>(&MAX_SIZE)->default_value(MAX_SIZE))
        ("help,h", po::value<bool>()->default_value(false)->implicit_value(true), "produce help")
        ("verbose,v", po::value<bool>(&VERBOSE)
            ->default_value(VERBOSE)->implicit_value(true), "verbose output")
        ("test", po::value<int>()->default_value(0)->implicit_value(1))
        ("mode", po::value<int>(&MODE)->default_value(MODE))
        ("out_dir,o", po::value<std::string>()->default_value(OUT_DIR.string()))
        ("seed", po::value<unsigned int>(&SEED)->default_value(SEED))
        ("nsam", po::value<size_t>(&NSAM)->default_value(NSAM))
        ("howmany", po::value<size_t>(&HOWMANY)->default_value(HOWMANY))
    ;
    return description;
}

//! Unit test for each class
inline void test() {HERE;
    Gland::unit_test();
    Tissue::unit_test();
}

Simulation::Simulation(int argc, char* argv[]) {HERE;
    std::vector<std::string> arguments(argv, argv + argc);
    COMMAND_ARGS = wtl::str_join(arguments, " ") ;

    std::ostringstream oss;
    oss << wtl::strftime("tumopp_%Y%m%d_%H%M_") << ::getpid();
    OUT_DIR = oss.str();

    namespace po = boost::program_options;
    po::options_description description;
    description.add(opt_description());
    description.add(Gland::opt_description());
    description.add(Tissue::opt_description());

    po::positional_options_description positional;
    positional.add("nsam", 1);
    positional.add("howmany", 1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
              options(description).
              positional(positional).run(), vm);
    po::notify(vm);

    if (vm["help"].as<bool>()) {
        description.print(std::cout);
        exit(0);
    }
    wtl::sfmt().seed(SEED);
    CONFIG_STRING = wtl::flags_into_string(description, vm);
    if (VERBOSE) {
        std::cerr << wtl::iso8601datetime() << std::endl;
        std::cerr << CONFIG_STRING << std::endl;
    }
    switch (vm["test"].as<int>()) {
      case 0:
        break;
      case 1:
        test();
        exit(0);
      default:
        exit(1);
    }
    OUT_DIR = fs::path(vm["out_dir"].as<std::string>());
    OUT_DIR = fs::system_complete(OUT_DIR);
}

void Simulation::run() const {HERE;
    std::cout << COMMAND_ARGS << "\n" << SEED << "\n";
    static size_t i = 0;
    Tissue tissue;

    const double expected_cs = tissue.coord_func()->cross_section(MAX_SIZE);
    derr("predicted section size: " << expected_cs << std::endl);
    if (NSAM > expected_cs) {
        std::cout.precision(1);
        std::cout << "ERROR: NSAM=" << NSAM
            << " is larger than the expected cross section size "
            << std::fixed << expected_cs << std::endl;
        exit(1);
    }

    switch (Tissue::DIMENSIONS()) {
      case 2: {
        tissue.grow(expected_cs);
        for (size_t i=0; i<HOWMANY; ++i) {
            tissue.write_segsites(std::cout, tissue.sample_random(NSAM));
        }
        break;
      }
      case 3: {
        tissue.grow(MAX_SIZE);
        for (size_t i=0; i<HOWMANY; ++i) {
            auto section = tissue.sample_if([](const std::vector<int>& coord) {
                return coord[2] == 0;  // z-axis
            });
            section = wtl::sample(section, NSAM, wtl::prandom());
            tissue.write_segsites(std::cout, section) << std::flush;
        }
        break;
      }
      default:
        exit(1);
    }

    if (VERBOSE) {
        auto SUB_DIR = OUT_DIR;
        SUB_DIR += wtl::strprintf("-%d", i);
        derr("mkdir && cd to " << SUB_DIR << std::endl);
        fs::create_directory(SUB_DIR);
        wtl::cd(SUB_DIR.string());
        wtl::Fout{"program_options.conf"} << CONFIG_STRING;
        wtl::gzip{wtl::Fout{"mutation_history.tsv.gz"}} << tissue.mutation_history();
        wtl::gzip{wtl::Fout{"population.tsv.gz"}}
            << tissue.snapshot_header() << tissue.snapshot();
        wtl::gzip{wtl::Fout{"evolution_history.tsv.gz"}}
            << tissue.snapshot_header() << tissue.evolution_history();
        std::cerr << wtl::iso8601datetime() << std::endl;
    }
    ++i;
}
