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
#include <wtl/exception.hpp>
#include <clippson/clippson.hpp>

namespace tumopp {

//! Global variables mapper of commane-line arguments
nlohmann::json VM;

//! Options description for general purpose
inline clipp::group general_options(nlohmann::json* vm) {
    return (
      wtl::option(vm, {"h", "help"}, false, "print this help"),
      wtl::option(vm, {"version"}, false, "print version"),
      wtl::option(vm, {"v", "verbose"}, false, "verbose output")
    ).doc("General:");
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
    `-o,--outdir`       | -              | -
    `-I,--interval`     | -              | -
    `-R,--record`       | -              | -
    `--seed`            | -              | -
*/
inline clipp::group simulation_options(nlohmann::json* vm) {
    const std::string OUT_DIR = wtl::strftime("tumopp_%Y%m%d_%H%M%S");
    const int seed = std::random_device{}(); // 32-bit signed integer for R
    return (
      wtl::option(vm, {"D", "dimensions"}, 3u),
      wtl::option(vm, {"C", "coord"}, "moore",
        "Coordinate/neighborhood system {neumann, moore, hex}"),
      wtl::option(vm, {"L", "local"}, "const",
        "E2 {const, step, linear}"),
      wtl::option(vm, {"P", "path"}, "random",
        "Push method {1: random, 2: roulette, 3: mindrag, 4: minstraight, 5: stroll}"),
      wtl::option(vm, {"O", "origin"}, 1u),
      wtl::option(vm, {"N", "max"}, 16384u),
      wtl::option(vm, {"T", "plateau"}, 0.0),
      wtl::option(vm, {"U", "mutate"}, 0u),
      wtl::option(vm, {"treatment"}, 0.0),
      wtl::option(vm, {"resistant"}, 3u),
      wtl::option(vm, {"o", "outdir"}, OUT_DIR),
      wtl::option(vm, {"I", "interval"}, 0.0),
      wtl::option(vm, {"R", "record"}, 0u),
      wtl::option(vm, {"extinction"}, 100u),
      wtl::option(vm, {"benchmark"}, false),
      wtl::option(vm, {"seed"}, seed)
    ).doc("Simulation:");
}

//! Parameters of Cell class
/*! @ingroup params

    Command line option | Symbol              | Variable                  |
    ------------------- | ------------------- | ------------------------- |
    `-b,--beta0`        | \f$\beta_0\f$       | EventRates::birth_rate
    `-d,--delta0`       | \f$\delta_0\f$      | EventRates::death_rate
    `-a,--alpha0`       | \f$\alpha_0\f$      | EventRates::death_prob
    `-m,--rho0`         | \f$\rho_0\f$        | EventRates::migra_rate
    `-k,--shape`        | \f$k\f$             | CellParams::GAMMA_SHAPE
    `-p,--symmetric`    | \f$p_s\f$           | CellParams::PROB_SYMMETRIC_DIVISION
    `-r,--prolif`       | \f$\omega_{\max}\f$ | CellParams::MAX_PROLIFERATION_CAPACITY
    `--ub`              | \f$\mu_\beta\f$     | CellParams::RATE_BIRTH
    `--ud`              | \f$\mu_\delta\f$    | CellParams::RATE_DEATH
    `--um`              | \f$\mu_\rho\f$      | CellParams::RATE_MIGRA
    `--mb`              | \f$\bar s_\beta\f$  | CellParams::MEAN_BIRTH
    `--md`              | \f$\bar s_\delta\f$ | CellParams::MEAN_DEATH
    `--mm`              | \f$\bar s_\rho\f$   | CellParams::MEAN_MIGRA
    `--sb`              | \f$\sigma_\beta\f$  | CellParams::SD_BIRTH
    `--sd`              | \f$\sigma_\delta\f$ | CellParams::SD_DEATH
    `--sm`              | \f$\sigma_\rho\f$   | CellParams::SD_MIGRA
*/
inline clipp::group
cell_options(nlohmann::json* vm, EventRates* init_event_rates, CellParams* cell_params) {
    return (
      wtl::option(vm, {"b", "beta0"}, &init_event_rates->birth_rate),
      wtl::option(vm, {"d", "delta0"}, &init_event_rates->death_rate),
      wtl::option(vm, {"a", "alpha0"}, &init_event_rates->death_prob),
      wtl::option(vm, {"m", "rho0"}, &init_event_rates->migra_rate),
      wtl::option(vm, {"k", "shape"}, &cell_params->GAMMA_SHAPE),
      wtl::option(vm, {"p", "symmetric"}, &cell_params->PROB_SYMMETRIC_DIVISION),
      wtl::option(vm, {"r", "prolif"}, &cell_params->MAX_PROLIFERATION_CAPACITY),
      wtl::option(vm, {"ub"}, &cell_params->RATE_BIRTH),
      wtl::option(vm, {"ud"}, &cell_params->RATE_DEATH),
      wtl::option(vm, {"ua"}, &cell_params->RATE_ALPHA),
      wtl::option(vm, {"um"}, &cell_params->RATE_MIGRA),
      wtl::option(vm, {"mb"}, &cell_params->MEAN_BIRTH),
      wtl::option(vm, {"md"}, &cell_params->MEAN_DEATH),
      wtl::option(vm, {"ma"}, &cell_params->MEAN_ALPHA),
      wtl::option(vm, {"mm"}, &cell_params->MEAN_MIGRA),
      wtl::option(vm, {"sb"}, &cell_params->SD_BIRTH),
      wtl::option(vm, {"sd"}, &cell_params->SD_DEATH),
      wtl::option(vm, {"sa"}, &cell_params->SD_ALPHA),
      wtl::option(vm, {"sm"}, &cell_params->SD_MIGRA)
    ).doc("Cell:");
}

Simulation::Simulation(const std::vector<std::string>& arguments)
: init_event_rates_(std::make_unique<EventRates>()),
  cell_params_(std::make_unique<CellParams>()) {
    std::ios::sync_with_stdio(false);
    std::cin.tie(0);
    std::cout.precision(9);
    std::cerr.precision(6);

    VM.clear();
    nlohmann::json vm_local;
    auto cli = (
      general_options(&vm_local),
      simulation_options(&VM),
      cell_options(&VM, init_event_rates_.get(), cell_params_.get())
    );
    wtl::parse(cli, arguments);
    auto fmt = wtl::doc_format();
    if (vm_local["help"]) {
        std::cout << "Usage: " << PROJECT_NAME << " [options]\n\n";
        std::cout << clipp::documentation(cli, fmt) << "\n";
        throw wtl::ExitSuccess();
    }
    if (vm_local["version"]) {
        std::cout << PROJECT_VERSION << "\n";
        throw wtl::ExitSuccess();
    }
    Cell::param(*cell_params_);
    config_ = VM.dump(2) + "\n";
    if (vm_local["verbose"]) {
        std::cerr << wtl::iso8601datetime() << std::endl;
        std::cerr << config_ << std::endl;
    }
}

Simulation::~Simulation() = default;

void Simulation::run() {
    const auto max_size = VM.at("max").get<size_t>();
    const double max_time = std::log2(max_size) * 100.0;
    const auto plateau_time = VM.at("plateau").get<double>();
    const auto treatment = VM.at("treatment").get<double>();
    const auto resistant = VM.at("resistant").get<size_t>();
    const auto allowed_extinction = VM.at("extinction").get<unsigned>();
    for (size_t i=0; i<allowed_extinction; ++i) {
        tissue_ = std::make_unique<Tissue>(
            VM.at("origin").get<size_t>(),
            VM.at("dimensions").get<unsigned>(),
            VM.at("coord").get<std::string>(),
            VM.at("local").get<std::string>(),
            VM.at("path").get<std::string>(),
            *init_event_rates_,
            VM.at("seed").get<uint_fast32_t>(),
            VM.at("benchmark").get<bool>()
        );
        bool success = tissue_->grow(
            max_size, max_time,
            VM.at("interval").get<double>(),
            VM.at("record").get<size_t>(),
            VM.at("mutate").get<size_t>()
        );
        if (success) break;
        std::cerr << "Trial " << i  << ": size = " << tissue_->size() << std::endl;
    }
    if (tissue_->size() != max_size) {
        std::cerr << "Warning: size = " << tissue_->size() << std::endl;
    }
    if (plateau_time > 0.0) {
        tissue_->plateau(plateau_time);
    }
    if (treatment > 0.0) {
        tissue_->treatment(treatment, resistant);
    }
}

std::string Simulation::outdir() const {
    return VM.at("outdir").get<std::string>();
}

std::string Simulation::history() const {
    std::ostringstream oss;
    tissue_->write_history(oss);
    return oss.str();
}
std::string Simulation::snapshots() const {
    if (!tissue_->has_snapshots()) return std::string{};
    std::ostringstream oss;
    tissue_->write_snapshots(oss);
    return oss.str();
}
std::string Simulation::drivers() const {
    if (!tissue_->has_drivers()) return std::string{};
    std::ostringstream oss;
    tissue_->write_drivers(oss);
    return oss.str();
}
std::string Simulation::benchmark() const {
    if (!tissue_->has_benchmark()) return std::string{};
    std::ostringstream oss;
    tissue_->write_benchmark(oss);
    return oss.str();
}

std::streambuf* std_cout_rdbuf(std::streambuf* buf) {
    return std::cout.rdbuf(buf);
}
std::streambuf* std_cerr_rdbuf(std::streambuf* buf) {
    return std::cerr.rdbuf(buf);
}

} // namespace tumopp
