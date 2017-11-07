/*! @file simulation.cpp
    @brief Implementation of Simulation class
    @defgroup params Parameters
*/
#include "simulation.hpp"
#include "cell.hpp"

#include <wtl/debug.hpp>
#include <wtl/iostr.hpp>
#include <wtl/getopt.hpp>
#include <wtl/zfstream.hpp>
#include <sfmt.hpp>

#include <boost/filesystem.hpp>

#include <cstdlib>
#include <algorithm>

namespace tumopp {

namespace fs = boost::filesystem;
namespace po = boost::program_options;

//! Options description for general purpose
inline po::options_description general_desc() {HERE;
    po::options_description description("General");
    description.add_options()
        ("help,h", po::bool_switch(), "print this help")
        ("verbose,v", po::bool_switch(), "verbose output")
        ("test", po::value<int>()->default_value(0)->implicit_value(1));
    return description;
}

//! Program options
/*! @ingroup params
    @return Program options description

    Command line option | Symbol         | Variable
    ------------------- | -------------- | -------------------------
    `-N,--max`          | \f$N_\max\f$   | Simulation::MAX_SIZE
    `-T,--plateau`      | -              | Simulation::PLATEAU_TIME
    `-w,--write`        | -              | Simulation::WRITE_TO_FILES
    `-o,--out_dir`      | -              | Simulation::OUT_DIR
    `--seed`            | -              | Simulation::SEED
*/
po::options_description Simulation::options_desc() {HERE;
    po::options_description description("Simulation");
    description.add_options()
        ("max,N", po::value(&MAX_SIZE)->default_value(MAX_SIZE))
        ("plateau,T", po::value(&PLATEAU_TIME)->default_value(PLATEAU_TIME))
        ("write,w", po::bool_switch(&WRITE_TO_FILES))
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

[[noreturn]] void Simulation::help_and_exit() {HERE;
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
        Cell::test();
        Tissue::test();
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
        if (tissue_.grow(MAX_SIZE, PLATEAU_TIME)) break;
    }
    if (tissue_.size() != MAX_SIZE) {
        std::cerr << "Warning: tissue_.size() " << tissue_.size() << std::endl;
    }
}

void Simulation::write() const {HERE;
    std::cout << "tumopp " << COMMAND_ARGS << "\n" << SEED << "\n";

    if (Tissue::DIMENSIONS() == 3U) {
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
        fs::current_path(OUT_DIR);
        wtl::opfstream{"program_options.conf"} << CONFIG_STRING;
        wtl::ozfstream{"population.tsv.gz"}
            << tissue_.specimens();
        wtl::ozfstream{"snapshots.tsv.gz"}
            << tissue_.snapshots();
        wtl::ozfstream{"drivers.tsv.gz"}
            << tissue_.drivers();
        wtl::ozfstream{"distance.tsv.gz"}
            << tissue_.pairwise_distance(std::min(200UL, tissue_.size() / 2U));
        std::cerr << wtl::iso8601datetime() << std::endl;
    }
}

} // namespace tumopp
