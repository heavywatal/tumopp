/*! @file simulation.cpp
    @brief Implementation of Simulation class
    @defgroup params Parameters
*/
#include "simulation.hpp"
#include "tissue.hpp"
#include "cell.hpp"
#include "version.hpp"

#include <wtl/debug.hpp>
#include <wtl/iostr.hpp>
#include <wtl/chrono.hpp>
#include <wtl/getopt.hpp>
#include <wtl/exception.hpp>
#include <wtl/zlib.hpp>
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
        ("version", po::bool_switch(), "print version")
        ("verbose,v", po::bool_switch(), "verbose output")
    ;
    return description;
}

//! Parameters of Simulation class
/*! @ingroup params

    Command line option | Symbol         | Variable                  |
    ------------------- | -------------- | ------------------------- |
    `-N,--max`          | \f$N_\max\f$   | -
    `-T,--plateau`      | -              | -
    `--npair`           | -              | -
    `-o,--outdir`       | -              | -
    `--seed`            | -              | -
*/
po::options_description Simulation::options_desc() {HERE;
    const std::string OUT_DIR = wtl::strftime("tumopp_%Y%m%d_%H%M%S");
    po::options_description description("Simulation");
    description.add_options()
      ("max,N", po::value<size_t>()->default_value(16384u))
      ("plateau,T", po::value<double>()->default_value(0.0))
      ("treatment", po::value<double>()->default_value(0.0))
      ("npair", po::value<size_t>()->default_value(0u))
      ("outdir,o", po::value<std::string>()->default_value(OUT_DIR))
      ("extinction", po::value<unsigned>()->default_value(100u))
      ("seed", po::value<uint32_t>()->default_value(std::random_device{}()));
    description.add(Cell::opt_description());
    description.add(Tissue::opt_description());
    return description;
}

//! Positional arguments to produce ms-like output
/*! @ingroup params

    Command line option | Symbol         | Variable                  |
    ------------------- | -------------- | ------------------------- |
    `nsam`              | -              | -
    `howmany`           | -              | -
*/
po::options_description Simulation::positional_desc() {HERE;
    po::options_description description("Positional");
    description.add_options()
        ("nsam", po::value<size_t>()->default_value(0u))
        ("howmany", po::value<size_t>()->default_value(1u));
    return description;
}

[[noreturn]] void Simulation::help_and_exit() {HERE;
    auto description = general_desc();
    description.add(options_desc());
    // do not print positional arguments as options
    std::cout << "Usage: tumopp [options] nsam howmany\n\n";
    description.print(std::cout);
    throw wtl::ExitSuccess();
}

[[noreturn]] void Simulation::version_and_exit() {HERE;
    std::cout << GIT_COMMIT_HASH
              << " [" << GIT_BRANCH << "] "
              << GIT_COMMIT_TIME << std::endl;
    throw wtl::ExitSuccess();
}

Simulation::Simulation(const std::vector<std::string>& arguments)
: vars_(std::make_unique<po::variables_map>()) {HERE;
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
    if (vm["version"].as<bool>()) {version_and_exit();}
    po::notify(vm);
    Cell::init_distributions();
    wtl::sfmt64().seed(vm["seed"].as<uint32_t>());

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
    const auto treatment = vm["treatment"].as<double>();
    const auto allowed_extinction = vm["extinction"].as<unsigned>();
    const double max_time = std::log2(max_size) * 100.0;

    for (size_t i=0; i<allowed_extinction; ++i) {
        tissue_ = std::make_unique<Tissue>();
        if (tissue_->grow(max_size, max_time)) break;
    }
    if (tissue_->size() != max_size) {
        std::cerr << "Warning: tissue_.size() " << tissue_->size() << std::endl;
    }
    if (plateau_time > 0.0) {
        tissue_->plateau(plateau_time);
    }
    if (treatment > 0.0) {
        tissue_->treatment(treatment);
    }
}

void Simulation::write() const {HERE;
    auto& vm = *vars_;
    const auto outdir = vm["outdir"].as<std::string>();
    const auto npair = vm["npair"].as<size_t>();
    if (outdir.empty()) return;
    DCERR("mkdir && cd to " << outdir << std::endl);
    fs::create_directory(outdir);
    fs::current_path(outdir);
    std::cerr << "Output: " << outdir << "\n";
    wtl::make_ofs("program_options.conf") << config_string_;
    wtl::zlib::ofstream{"population.tsv.gz"} << tissue_->specimens();
    if (!tissue_->snapshots_empty()) {
        wtl::zlib::ofstream{"snapshots.tsv.gz"} << tissue_->snapshots();
    }
    if (!tissue_->drivers_empty()) {
        wtl::zlib::ofstream{"drivers.tsv.gz"} << tissue_->drivers();
    }
    if (npair > 0u) {
        wtl::zlib::ofstream{"distance.tsv.gz"} << tissue_->pairwise_distance(npair);
    }
}

void Simulation::ms(std::ostream& ost) const {HERE;
    auto& vm = *vars_;
    const auto nsam = vm["nsam"].as<size_t>();
    const auto howmany = vm["howmany"].as<size_t>();
    const auto seed = vm["seed"].as<uint32_t>();
    if ((nsam < 1u) || (howmany < 1u)) return;
    ost << "tumopp " << command_args_ << "\n" << seed << "\n";

    if (Tissue::DIMENSIONS() == 3U) {
        for (size_t i=0; i<howmany; ++i) {
            tissue_->write_segsites(ost, tissue_->sample_section(nsam));
        }
    } else {
        for (size_t i=0; i<howmany; ++i) {
            tissue_->write_segsites(ost, tissue_->sample_random(nsam));
        }
    }
}

std::string Simulation::specimens() const {return tissue_->specimens();}
std::string Simulation::snapshots() const {return tissue_->snapshots();}
std::string Simulation::drivers() const {return tissue_->drivers();}
std::string Simulation::pairwise_distance(size_t npair) const {return tissue_->pairwise_distance(npair);}
std::string Simulation::ms() const {
    std::ostringstream oss;
    ms(oss);
    return oss.str();
}

} // namespace tumopp
