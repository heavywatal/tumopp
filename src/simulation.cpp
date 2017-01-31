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

namespace tumopp {

namespace fs = boost::filesystem;
namespace po = boost::program_options;

inline po::options_description general_desc() {HERE;
    po::options_description description("General");
    description.add_options()
        ("help,h", po::value<bool>()->default_value(false)->implicit_value(true), "print this help")
        ("verbose,v", po::value<bool>()->default_value(false)->implicit_value(true), "verbose output")
        ("test", po::value<int>()->default_value(0)->implicit_value(1));
    return description;
}

//! Program options
/*! @return Program options description

    Command line option | Symbol         | Variable
    --------------------| -------------- | -------------------------
    `-N,--max`          | -              | Simulation::MAX_SIZE
*/
po::options_description Simulation::options_desc() {HERE;
    po::options_description description("Simulation");
    description.add_options()
        ("max,N", po::value(&MAX_SIZE)->default_value(MAX_SIZE))
        ("write,w", po::value(&WRITE_TO_FILES)->default_value(WRITE_TO_FILES)->implicit_value(true))
        ("out_dir,o", po::value(&OUT_DIR)->default_value(OUT_DIR))
        ("seed", po::value(&SEED)->default_value(SEED));
    description.add(Cell::opt_description());
    description.add(Tissue::opt_description());
    return description;
}

po::options_description Simulation::positional_desc() {HERE;
    po::options_description description("Positional");
    description.add_options()
        ("nsam", po::value(&NSAM)->default_value(NSAM))
        ("howmany", po::value(&HOWMANY)->default_value(HOWMANY));
    return description;
}

void Simulation::help_and_exit() {HERE;
    auto description = general_desc();
    description.add(options_desc());
    // do not print positional arguments as options
    std::cout << "Usage: tumopp [options] nsam howmany\n" << std::endl;
    description.print(std::cout);
    throw wtl::ExitSuccess();
}

//! Unit test for each class
inline void test(const int flg) {HERE;
    switch (flg) {
      case 0:
        break;
      case 1:
        Cell::unit_test();
        Tissue::unit_test();
        throw wtl::ExitSuccess();
      default:
        throw std::runtime_error("Unknown argument for --test");
    }
}

Simulation::Simulation(const std::vector<std::string>& arguments) {HERE;
    std::ios::sync_with_stdio(false);
    std::cin.tie(0);
    std::cout.precision(15);
    std::cerr.precision(6);
    COMMAND_ARGS = wtl::str_join(arguments, " ");
    OUT_DIR = wtl::strftime("tumopp_%Y%m%d_%H%M_") + std::to_string(::getpid());

    auto description = general_desc();
    description.add(options_desc());
    description.add(positional_desc());
    po::positional_options_description positional;
    positional.add("nsam", 1)
              .add("howmany", 1);
    po::variables_map vm;
    po::store(po::command_line_parser(arguments).
              options(description).
              positional(positional).run(), vm);
    if (vm["help"].as<bool>()) {help_and_exit();}
    po::notify(vm);
    Cell::init_distributions();
    wtl::sfmt().seed(SEED);

    CONFIG_STRING = wtl::flags_into_string(vm);
    if (vm["verbose"].as<bool>()) {
        std::cerr << wtl::iso8601datetime() << std::endl;
        std::cerr << CONFIG_STRING << std::endl;
    }
    test(vm["test"].as<int>());
    if (NSAM > vm["max"].as<size_t>()) {
        std::ostringstream oss;
        oss << "NSAM=" << NSAM
            << " is larger than final tumor size "
            << vm["max"].as<size_t>();
        throw std::runtime_error(oss.str());
    }
    OUT_DIR = fs::system_complete(OUT_DIR).string();
}

void Simulation::run() {HERE;
    for (size_t i=0; i<10; ++i) {
        if (tissue_.grow(MAX_SIZE)) break;
    }
    if (tissue_.size() != MAX_SIZE) {
        std::cerr << "Warning: tissue_.size() " << tissue_.size() << std::endl;
    }
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

    if (WRITE_TO_FILES) {
        DCERR("mkdir && cd to " << OUT_DIR << std::endl);
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
