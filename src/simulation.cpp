/*! @file simulation.cpp
    @brief Implementation of Simulation class
    @defgroup params Parameters
*/
#include "simulation.hpp"
#include "tissue.hpp"
#include "cell.hpp"

#include <wtl/debug.hpp>
#include <wtl/iostr.hpp>
#include <wtl/chrono.hpp>
#include <wtl/getopt.hpp>
#include <wtl/exception.hpp>
#include <wtl/zfstream.hpp>
#include <wtl/random.hpp>
#include <sfmt.hpp>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

namespace tumopp {

namespace fs = boost::filesystem;
namespace po = boost::program_options;

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
    `-N,--max`          | \f$N_\max\f$   |
    `-T,--plateau`      | -              |
    `-o,--outdir`       | -              |
    `--seed`            | -              |
*/
po::options_description Simulation::options_desc() {HERE;
    const std::string OUT_DIR = wtl::strftime("tumopp_%Y%m%d_%H%M%S");
    po::options_description description("Simulation");
    description.add_options()
      ("max,N", po::value<size_t>()->default_value(16384u))
      ("plateau,T", po::value<double>()->default_value(0.0))
      ("outdir,o", po::value<std::string>()->default_value("")->implicit_value(OUT_DIR))
      ("seed", po::value<size_t>()->default_value(wtl::random_device_64{}()));
    description.add(Cell::opt_description());
    description.add(Tissue::opt_description());
    return description;
}

po::options_description Simulation::positional_desc() {HERE;
    po::options_description description("Positional");
    description.add_options()
        ("nsam", po::value<size_t>()->default_value(20u))
        ("howmany", po::value<size_t>()->default_value(1u));
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
: tissue_(std::make_unique<Tissue>()),
  vars_(std::make_unique<po::variables_map>()) {HERE;
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
    auto& vm = *vars_;
    po::store(po::command_line_parser(arguments).
              options(description).
              positional(positional).run(), vm);
    if (vm["help"].as<bool>()) {help_and_exit();}
    po::notify(vm);
    Cell::init_distributions();
    wtl::sfmt64().seed(vm["seed"].as<size_t>());

    config_string_ = wtl::flags_into_string(vm);
    if (vm["verbose"].as<bool>()) {
        std::cerr << wtl::iso8601datetime() << std::endl;
        std::cerr << config_string_ << std::endl;
    }
    if (vm["nsam"].as<size_t>() > vm["max"].as<size_t>()) {
        std::ostringstream oss;
        oss << "nsam_=" << vm["nsam"].as<size_t>()
            << " is larger than final tumor size "
            << vm["max"].as<size_t>();
        throw std::runtime_error(oss.str());
    }
}

Simulation::~Simulation() = default;

void Simulation::run() {HERE;
    auto& vm = *vars_;
    const auto max_size = vm["max"].as<size_t>();
    const auto plateau_time = vm["plateau"].as<double>();

    for (size_t i=0; i<10u; ++i) {
        if (tissue_->grow(max_size, plateau_time)) break;
    }
    if (tissue_->size() != max_size) {
        std::cerr << "Warning: tissue_.size() " << tissue_->size() << std::endl;
    }
}

void Simulation::write() const {HERE;
    auto& vm = *vars_;
    const auto nsam = vm["nsam"].as<size_t>();
    const auto howmany = vm["howmany"].as<size_t>();
    const auto outdir = vm["outdir"].as<std::string>();
    const auto seed = vm["seed"].as<size_t>();
    std::cout << "tumopp " << command_args_ << "\n" << seed << "\n";

    if (Tissue::DIMENSIONS() == 3U) {
        for (size_t i=0; i<howmany; ++i) {
            tissue_->write_segsites(std::cout, tissue_->sample_section(nsam));
        }
    } else {
        for (size_t i=0; i<howmany; ++i) {
            tissue_->write_segsites(std::cout, tissue_->sample_random(nsam));
        }
    }

    if (!outdir.empty()) {
        DCERR("mkdir && cd to " << outdir << std::endl);
        fs::create_directory(outdir);
        fs::current_path(outdir);
        wtl::make_ofs("program_options.conf") << config_string_;
        wtl::ozfstream{"population.tsv.gz"}
            << tissue_->specimens();
        wtl::ozfstream{"snapshots.tsv.gz"}
            << tissue_->snapshots();
        wtl::ozfstream{"drivers.tsv.gz"}
            << tissue_->drivers();
        wtl::ozfstream{"distance.tsv.gz"}
            << tissue_->pairwise_distance(std::min(200UL, tissue_->size() / 2U));
        std::cerr << wtl::iso8601datetime() << std::endl;
    }
}

std::vector<std::string> Simulation::results(size_t npair) const {
    return {
        config_string_,
        tissue_->specimens(),
        tissue_->drivers(),
        tissue_->pairwise_distance(npair)
    };
}

} // namespace tumopp
