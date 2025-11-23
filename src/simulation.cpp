/*! @file simulation.cpp
    @brief Implementation of Simulation class
    @defgroup params Parameters
*/
#include "simulation.hpp"
#include "tissue.hpp"
#include "cell.hpp"
#include "random.hpp"
#include "version.hpp"

#include <fmt/format.h>
#include <wtl/zlib.hpp>
#include <wtl/chrono.hpp>
#include <wtl/signed.hpp>
#include <clippson/clippson.hpp>

#include <filesystem>
#include <fstream>
#include <cstdio>

namespace tumopp {

//! Global variables mapper of command-line arguments
nlohmann::json VM;

//! Options description for general purpose
inline clipp::group general_options(nlohmann::json* vm) {
    return (
      clippson::option(vm, {"h", "help"}, false, "Print this help"),
      clippson::option(vm, {"version"}, false, "Print version")
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
    return (
      clippson::option(vm, {"D", "dimensions"}, 3),
      clippson::option(vm, {"C", "coord"},
        "moore",
        "Neighborhood"
        " {neumann, moore, hex}"), // TODO
      clippson::option(vm, {"L", "local"},
        "const",
        "E2: resource competition"
        " {const, step, linear}"), // TODO
      clippson::option(vm, {"P", "path"},
        "random",
        "Push method"
        " {random, roulette, mindrag, minstraight, stroll}"), // TODO
      clippson::option(vm, {"O", "origin"}, 1),
      clippson::option(vm, {"N", "max"}, 16384,
        "Maximum number of cells to simulate"),
      clippson::option(vm, {"max_time"}, 0.0,
        "Maximum time to simulate"),
      clippson::option(vm, {"T", "plateau"}, 0.0,
        "Duration of turn-over phase after population growth"),
      clippson::option(vm, {"U", "mutate"}, 0,
        "Introduce a driver mutation to U-th cell"),
      clippson::option(vm, {"treatment"}, 0.0),
      clippson::option(vm, {"resistant"}, 3),
      clippson::option(vm, {"o", "outdir"}, "tumopp_%Y%m%d_%H%M%S"),
      clippson::option(vm, {"I", "interval"}, 0.0,
        "Time interval to take snapshots"),
      clippson::option(vm, {"R", "record"}, 0,
        "Tumor size to stop taking snapshots"),
      clippson::option(vm, {"extinction"}, 100,
        "Maximum number of trials in case of extinction"),
      clippson::option(vm, {"benchmark"}, false),
      clippson::option(vm, {"seed"}, 0,
        "Initialized with `std::random_device` if 0"),
      clippson::option(vm, {"v", "verbose"}, false, "Verbose output")
    ).doc("Simulation:");
}

//! Parameters of Cell class
/*! @ingroup params

    Command line option | Symbol              | Variable                  |
    ------------------- | ------------------- | ------------------------- |
    `-b,--beta0`        | \f$\beta_0\f$       | EventRates::birth_rate
    `-d,--delta0`       | \f$\delta_0\f$      | EventRates::death_rate
    `-a,--alpha0`       | \f$\alpha_0\f$      | EventRates::death_prob
    `-m,--rho0`         | \f$\rho_0\f$        | EventRates::migration_rate
    `-k,--shape`        | \f$k\f$             | CellParams::GAMMA_SHAPE
    `-p,--symmetric`    | \f$p_s\f$           | CellParams::PROB_SYMMETRIC_DIVISION
    `-r,--prolif`       | \f$\omega_{\max}\f$ | CellParams::MAX_PROLIFERATION_CAPACITY
    `--ub`              | \f$\mu_\beta\f$     | CellParams::RATE_BIRTH
    `--ud`              | \f$\mu_\delta\f$    | CellParams::RATE_DEATH
    `--um`              | \f$\mu_\rho\f$      | CellParams::RATE_MIG
    `--mb`              | \f$\bar s_\beta\f$  | CellParams::MEAN_BIRTH
    `--md`              | \f$\bar s_\delta\f$ | CellParams::MEAN_DEATH
    `--mm`              | \f$\bar s_\rho\f$   | CellParams::MEAN_MIG
    `--sb`              | \f$\sigma_\beta\f$  | CellParams::SD_BIRTH
    `--sd`              | \f$\sigma_\delta\f$ | CellParams::SD_DEATH
    `--sm`              | \f$\sigma_\rho\f$   | CellParams::SD_MIG
*/
inline clipp::group
cell_options(nlohmann::json* vm, EventRates* init_event_rates, CellParams* cell_params) {
    return (
      clippson::option(vm, {"b", "beta0"}, &init_event_rates->birth_rate, "Basic birth rate"),
      clippson::option(vm, {"d", "delta0"}, &init_event_rates->death_rate, "Basic death rate"),
      clippson::option(vm, {"a", "alpha0"}, &init_event_rates->death_prob,
        "Basic death rate on cell division attempt"),
      clippson::option(vm, {"m", "rho0"}, &init_event_rates->migration_rate, "Basic migration rate"),
      clippson::option(vm, {"k", "shape"}, &cell_params->GAMMA_SHAPE,
        "Shape parameter of waiting time distribution for cell division"),
      clippson::option(vm, {"p", "symmetric"}, &cell_params->PROB_SYMMETRIC_DIVISION,
        "p_s: Probability of symmetric division"),
      clippson::option(vm, {"r", "prolif"}, &cell_params->MAX_PROLIFERATION_CAPACITY,
        "ω: Maximum number of division for a TAC"),
      (
        clippson::option(vm, {"ub"}, &cell_params->RATE_BIRTH, "μ_β"),
        clippson::option(vm, {"ud"}, &cell_params->RATE_DEATH, "μ_δ"),
        clippson::option(vm, {"ua"}, &cell_params->RATE_ALPHA, "μ_α"),
        clippson::option(vm, {"um"}, &cell_params->RATE_MIG, "μ_ρ")
      ).doc("Rate of driver mutations:"),
      (
        clippson::option(vm, {"mb"}, &cell_params->MEAN_BIRTH, "E[s_β]"),
        clippson::option(vm, {"md"}, &cell_params->MEAN_DEATH, "E[s_δ]"),
        clippson::option(vm, {"ma"}, &cell_params->MEAN_ALPHA, "E[s_α]"),
        clippson::option(vm, {"mm"}, &cell_params->MEAN_MIG, "E[s_ρ]")
      ).doc("Mean effect of driver mutations:"),
      (
        clippson::option(vm, {"sb"}, &cell_params->SD_BIRTH, "σ_β"),
        clippson::option(vm, {"sd"}, &cell_params->SD_DEATH, "σ_δ"),
        clippson::option(vm, {"sa"}, &cell_params->SD_ALPHA, "σ_α"),
        clippson::option(vm, {"sm"}, &cell_params->SD_MIG, "σ_ρ")
      ).doc("SD of driver mutations:")
    ).doc("Cell:");
}

Simulation::Simulation(const std::vector<std::string>& arguments)
: init_event_rates_(std::make_unique<EventRates>()),
  cell_params_(std::make_unique<CellParams>()) {
    std::ios::sync_with_stdio(false);
    std::cin.tie(0);
    VM.clear();
    nlohmann::json vm_local;
    auto cli = (
      general_options(&vm_local),
      simulation_options(&VM),
      cell_options(&VM, init_event_rates_.get(), cell_params_.get())
    );
    clippson::parse(cli, arguments);
    if (vm_local["help"]) {
        auto doc_fmt = clippson::doc_format().paragraph_spacing(0);
        fmt::println("Usage: {} [options]\n", PROJECT_NAME);
        fmt::println("{}", clipp::documentation(cli, doc_fmt).str());
        throw exit_success();
    }
    if (vm_local["version"]) {
        fmt::println("{}", PROJECT_VERSION);
        throw exit_success();
    }
    Cell::param(*cell_params_);
    const auto vm_o = VM.at("outdir").get<std::string>();
    if (vm_o.find('%') != std::string::npos) {
      VM["outdir"] = wtl::strftime(vm_o.c_str());
    }
    if (VM.at("seed").get<int32_t>() == 0) {// 32-bit signed integer for R
      VM["seed"] = static_cast<int32_t>(std::random_device{}());
    }
    config_ = VM.dump(2) + "\n";
}

Simulation::~Simulation() = default;

void Simulation::run() {
    const auto max_size = VM.at("max").get<ptrdiff_t>();
    const double max_time = VM.at("max_time").get<double>();
    const auto plateau_time = VM.at("plateau").get<double>();
    const auto treatment = VM.at("treatment").get<double>();
    const auto resistant = VM.at("resistant").get<ptrdiff_t>();
    const auto allowed_extinction = VM.at("extinction").get<int>();
    urbg_t seeder(VM.at("seed").get<uint32_t>());
    for (int i=0; i<allowed_extinction; ++i) {
        tissue_ = std::make_unique<Tissue>(
            VM.at("origin").get<ptrdiff_t>(),
            VM.at("dimensions").get<int>(),
            VM.at("coord").get<std::string>(),
            VM.at("local").get<std::string>(),
            VM.at("path").get<std::string>(),
            *init_event_rates_,
            seeder(),
            VM.at("verbose").get<bool>(),
            VM.at("benchmark").get<bool>()
        );
        bool success = tissue_->grow(
            max_size,
            max_time > 0.0 ? max_time : std::log2(max_size) * 100.0,
            VM.at("interval").get<double>(),
            VM.at("record").get<ptrdiff_t>(),
            VM.at("mutate").get<ptrdiff_t>()
        );
        if (success) break;
        fmt::println(stderr, "Trial {}: size = {}", i, wtl::ssize(*tissue_));
    }
    if (max_time == 0.0 && wtl::ssize(*tissue_) != max_size) {
        fmt::println(stderr, "Warning: size = {}", wtl::ssize(*tissue_));
    }
    if (max_time == 0.0 && plateau_time > 0.0) {
        tissue_->plateau(plateau_time);
    }
    if (max_time == 0.0 && treatment > 0.0) {
        const ptrdiff_t margin = 10 * resistant + 10;
        tissue_->treatment(treatment, resistant);
        tissue_->grow(
            wtl::ssize(*tissue_) + margin,
            std::numeric_limits<double>::max(),
            VM.at("interval").get<double>()
        );
    }
}

//! Write config and simulation result to files
void Simulation::write() const {
    namespace fs = std::filesystem;
    const auto& outdir = VM.at("outdir").get<fs::path>();
    if (outdir.empty()) return;
    fs::create_directory(outdir);
    std::ofstream{outdir / "config.json"} << config_;
    {
        wtl::zlib::ofstream ofs{outdir / "population.tsv.gz"};
        ofs.exceptions(std::ios_base::failbit | std::ios_base::badbit);
        ofs << tissue_->str_history();
    }
    if (tissue_->has_snapshots()) {
        wtl::zlib::ofstream ofs{outdir / "snapshots.tsv.gz"};
        ofs.exceptions(std::ios_base::failbit | std::ios_base::badbit);
        ofs << tissue_->str_snapshots();
    }
    if (tissue_->has_drivers()) {
        wtl::zlib::ofstream ofs{outdir / "drivers.tsv.gz"};
        ofs.exceptions(std::ios_base::failbit | std::ios_base::badbit);
        ofs << tissue_->str_drivers();
    }
    if (!tissue_->str_benchmark().empty()) {
        wtl::zlib::ofstream ofs{outdir / "benchmark.tsv.gz"};
        ofs.exceptions(std::ios_base::failbit | std::ios_base::badbit);
        tissue_->benchmark_append(wtl::ssize(*tissue_) + 1);
        ofs << tissue_->str_benchmark();
    }
}

} // namespace tumopp
