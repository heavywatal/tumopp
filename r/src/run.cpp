// [[Rcpp::plugins(cpp14)]]
#include <Rcpp.h>
#include <tumopp/simulation.hpp>

//' Run C++ simulation
//' @return conf and population as strings
//' @rdname tumopp
// [[Rcpp::export]]
std::vector<std::string> cpp_tumopp(const std::vector<std::string>& args, unsigned int npair=0, unsigned int nsam=0) {
    try {
        tumopp::Simulation simulation(args);
        simulation.run();
        std::vector<std::string> output;
        output.reserve(5u);
        output.emplace_back(simulation.config_string());
        output.emplace_back(simulation.specimens());
        output.emplace_back(simulation.drivers());
        output.emplace_back(simulation.pairwise_distance(npair));
        if (nsam > 0u) {
          output.emplace_back(simulation.ms());
        }
        return output;
    } catch (const std::runtime_error& e) {
        Rcpp::Rcerr << e.what() << "\n";
    }
    return {};
}
