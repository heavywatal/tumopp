/*! @file simulation.cpp
    @brief Implementation of Simulation class
    @defgroup params Parameters
*/
#include "simulation.hpp"
#include "cell.hpp"

#include <wtl/debug.hpp>
#include <wtl/iostr.hpp>
#include <wtl/chrono.hpp>
#include <wtl/getopt.hpp>
#include <wtl/exception.hpp>
#include <wtl/zfstream.hpp>
#include <sfmt.hpp>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

namespace tumopp {

namespace fs = boost::filesystem;
namespace po = boost::program_options;

//! predefined name of output directory for -w option
const std::string OUT_DIR = wtl::strftime("tumopp_%Y%m%d_%H%M_") + std::to_string(::getpid());

//! Options description for general purpose
inline po::options_description general_desc() {HERE;
    po::options_description description("General");
    description.add_options()
        ("help,h", po::bool_switch(), "print this help")
        ("verbose,v", po::bool_switch(), "verbose output")
    ;
    return description;
}

//! Program options
/*! @ingroup params
    @return Program options description

    Command line option | Symbol         | Variable
    ------------------- | -------------- | -------------------------
    `-N,--max`          | \f$N_\max\f$   | Simulation::max_size_
    `-T,--plateau`      | -              | Simulation::plateau_time_
    `-o,--outdir`       | -              | Simulation::out_dir_
    `--seed`            | -              | Simulation::seed_
*/
po::options_description Simulation::options_desc() {HERE;
    po::options_description description("Simulation");
    description.add_options()
      ("max,N", po::value(&max_size_)->default_value(max_size_))
      ("plateau,T", po::value(&plateau_time_)->default_value(plateau_time_))
      ("outdir-predefined,w", po::bool_switch(), OUT_DIR.c_str())
      ("outdir,o", po::value(&out_dir_))
      ("seed", po::value(&seed_)->default_value(seed_));
    description.add(Cell::opt_description());
    description.add(Tissue::opt_description());
    return description;
}

po::options_description Simulation::positional_desc() {HERE;
    po::options_description description("Positional");
    description.add_options()
        ("nsam", po::value(&nsam_)->default_value(nsam_))
        ("howmany", po::value(&howmany_)->default_value(howmany_));
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

Simulation::Simulation(const std::vector<std::string>& arguments)
: seed_(std::random_device{}()) {HERE;
    std::ios::sync_with_stdio(false);
    std::cin.tie(0);
    std::cout.precision(15);
    std::cerr.precision(6);
    command_args_ = wtl::str_join(arguments, " ");

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
    wtl::sfmt64().seed(seed_);

    config_string_ = wtl::flags_into_string(vm);
    if (vm["verbose"].as<bool>()) {
        std::cerr << wtl::iso8601datetime() << std::endl;
        std::cerr << config_string_ << std::endl;
    }
    if (nsam_ > vm["max"].as<size_t>()) {
        std::ostringstream oss;
        oss << "nsam_=" << nsam_
            << " is larger than final tumor size "
            << vm["max"].as<size_t>();
        throw std::runtime_error(oss.str());
    }
    if (vm["outdir-predefined"].as<bool>()) {
        if (vm.count("outdir")) {
            throw std::runtime_error("cannot use -o and -w at the same time");
        }
        out_dir_ = fs::system_complete(OUT_DIR).string();
    }
}

void Simulation::run() {HERE;
    for (size_t i=0; i<10; ++i) {
        if (tissue_.grow(max_size_, plateau_time_)) break;
    }
    if (tissue_.size() != max_size_) {
        std::cerr << "Warning: tissue_.size() " << tissue_.size() << std::endl;
    }
}

void Simulation::write() const {HERE;
    std::cout << "tumopp " << command_args_ << "\n" << seed_ << "\n";

    if (Tissue::DIMENSIONS() == 3U) {
        for (size_t i=0; i<howmany_; ++i) {
            tissue_.write_segsites(std::cout, tissue_.sample_section(nsam_));
        }
    } else {
        for (size_t i=0; i<howmany_; ++i) {
            tissue_.write_segsites(std::cout, tissue_.sample_random(nsam_));
        }
    }

    if (!out_dir_.empty()) {
        DCERR("mkdir && cd to " << out_dir_ << std::endl);
        fs::create_directory(out_dir_);
        fs::current_path(out_dir_);
        wtl::make_ofs("program_options.conf") << config_string_;
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
