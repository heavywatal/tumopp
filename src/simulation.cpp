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
    `-D,--dimensions`   | -              | -
    `-C,--coord`        | -              | -
    `-L,--local`        | \f$E_2\f$      | -
    `-P,--path`         | -              | -
    `-O,--origin`       | \f$N_0\f$      | -
    `-N,--max`          | \f$N_\max\f$   | -
    `-T,--plateau`      | -              | -
    `-U,--mutate`       | \f$N_\mu\f$    | -
    `--npair`           | -              | -
    `-o,--outdir`       | -              | -
    `-I,--interval`     | -              | -
    `-R,--record`       | -              | -
    `--seed`            | -              | -
*/
po::options_description Simulation::options_desc() {HERE;
    const std::string OUT_DIR = wtl::strftime("tumopp_%Y%m%d_%H%M%S");
    po::options_description description("Simulation");
    description.add_options()
      ("dimensions,D", po::value<unsigned>()->default_value(3u))
      ("coord,C", po::value<std::string>()->default_value("moore"),
       "Coordinate/neighborhood system {neumann, moore, hex}")
      ("local,L", po::value<std::string>()->default_value("const"),
       "E2 {const, step, linear}")
      ("path,P", po::value<std::string>()->default_value("random"),
       "Push method {1: random, 2: roulette, 3: mindrag, 4: minstraight, 5: stroll}")
      ("origin,O", po::value<size_t>()->default_value(1u))
      ("max,N", po::value<size_t>()->default_value(16384u))
      ("plateau,T", po::value<double>()->default_value(0.0))
      ("mutate,U", po::value<size_t>()->default_value(std::numeric_limits<size_t>::max()))
      ("treatment", po::value<double>()->default_value(0.0))
      ("resistant", po::value<size_t>()->default_value(3u))
      ("npair", po::value<size_t>()->default_value(0u))
      ("outdir,o", po::value<std::string>()->default_value(OUT_DIR))
      ("interval,I", po::value<double>()->default_value(std::numeric_limits<double>::infinity()))
      ("record,R", po::value<size_t>()->default_value(0u))
      ("extinction", po::value<unsigned>()->default_value(100u))
      ("seed", po::value<uint32_t>()->default_value(std::random_device{}()));
    description.add(Cell::opt_description());
    description.add(cell_options());
    return description;
}

//! Parameters of Cell class
/*! @ingroup params

    Command line option | Symbol              | Variable                  |
    ------------------- | ------------------- | ------------------------- |
    `-b,--beta0`        | \f$\beta_0\f$       | EventRates::birth_rate
    `-d,--delta0`       | \f$\delta_0\f$      | EventRates::death_rate
    `-a,--alpha0`       | \f$\alpha_0\f$      | EventRates::death_prob
    `-m,--rho0`         | \f$\rho_0\f$        | EventRates::migra_rate
    `--ub`              | \f$\mu_\beta\f$     | DriverParams::RATE_BIRTH
    `--ud`              | \f$\mu_\delta\f$    | DriverParams::RATE_DEATH
    `--um`              | \f$\mu_\rho\f$      | DriverParams::RATE_MIGRA
    `--mb`              | \f$\bar s_\beta\f$  | DriverParams::MEAN_BIRTH
    `--md`              | \f$\bar s_\delta\f$ | DriverParams::MEAN_DEATH
    `--mm`              | \f$\bar s_\rho\f$   | DriverParams::MEAN_MIGRA
    `--sb`              | \f$\sigma_\beta\f$  | DriverParams::SD_BIRTH
    `--sd`              | \f$\sigma_\delta\f$ | DriverParams::SD_DEATH
    `--sm`              | \f$\sigma_\rho\f$   | DriverParams::SD_MIGRA
*/
po::options_description Simulation::cell_options() {HERE;
    init_event_rates_ = std::make_unique<EventRates>();
    driver_params_ = std::make_unique<DriverParams>();
    namespace po = boost::program_options;
    po::options_description desc{"Cell"};
    auto po_value = [](auto* var) {return po::value(var)->default_value(*var);};
    desc.add_options()
      ("beta0,b", po_value(&init_event_rates_->birth_rate))
      ("delta0,d", po_value(&init_event_rates_->death_rate))
      ("alpha0,a", po_value(&init_event_rates_->death_prob))
      ("rho0,m", po_value(&init_event_rates_->migra_rate))
      ("ub", po_value(&driver_params_->RATE_BIRTH))
      ("ud", po_value(&driver_params_->RATE_DEATH))
      ("um", po_value(&driver_params_->RATE_MIGRA))
      ("mb", po_value(&driver_params_->MEAN_BIRTH))
      ("md", po_value(&driver_params_->MEAN_DEATH))
      ("mm", po_value(&driver_params_->MEAN_MIGRA))
      ("sb", po_value(&driver_params_->SD_BIRTH))
      ("sd", po_value(&driver_params_->SD_DEATH))
      ("sm", po_value(&driver_params_->SD_MIGRA))
    ;
    return desc;
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
    Cell::init_distributions(*driver_params_);
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
    const auto dimensions = vm["dimensions"].as<unsigned>();
    const auto coord = vm["coord"].as<std::string>();
    const auto local = vm["local"].as<std::string>();
    const auto path = vm["path"].as<std::string>();
    const auto init_size = vm["origin"].as<size_t>();
    const auto max_size = vm["max"].as<size_t>();
    const auto plateau_time = vm["plateau"].as<double>();
    const auto mutate = vm["mutate"].as<size_t>();
    const auto treatment = vm["treatment"].as<double>();
    const auto resistant = vm["resistant"].as<size_t>();
    const auto interval = vm["interval"].as<double>();
    const auto record = vm["record"].as<size_t>();
    const auto allowed_extinction = vm["extinction"].as<unsigned>();
    const double max_time = std::log2(max_size) * 100.0;

    for (size_t i=0; i<allowed_extinction; ++i) {
        tissue_ = std::make_unique<Tissue>(init_size, dimensions, coord, local, path, *init_event_rates_);
        if (tissue_->grow(max_size, max_time, interval, record, mutate)) break;
    }
    if (tissue_->size() != max_size) {
        std::cerr << "Warning: tissue_.size() " << tissue_->size() << std::endl;
    }
    if (plateau_time > 0.0) {
        tissue_->plateau(plateau_time);
    }
    if (treatment > 0.0) {
        tissue_->treatment(treatment, resistant);
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

    if (tissue_->dimensions() == 3U) {
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
std::string Simulation::pairwise_distance(size_t npair) const {
    if (npair == 0u) npair = (*vars_)["npair"].as<size_t>();
    return tissue_->pairwise_distance(npair);
}
std::string Simulation::ms() const {
    std::ostringstream oss;
    ms(oss);
    return oss.str();
}

} // namespace tumopp
