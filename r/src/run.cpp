// [[Rcpp::plugins(cpp14)]]
#include <Rcpp.h>
#include <tumopp/simulation.hpp>

//' Run C++ simulation
//' @return conf and population as strings
//' @rdname tumopp
// [[Rcpp::export]]
std::vector<std::string> cpp_tumopp(const std::vector<std::string>& args, size_t npair=0) {
    try {
        tumopp::Simulation simulation(args);
        simulation.run();
        return simulation.results(npair);
    } catch (const std::runtime_error& e) {
        Rcpp::Rcerr << e.what() << "\n";
    }
    return {};
}
