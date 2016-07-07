// [[Rcpp::plugins(cpp14)]]
#include <Rcpp.h>
#include <tumopp/simulation.hpp>

//' Run C++ simulation
//' @param args string vector
//' @return conf and population as strings
//' @rdname cpp
// [[Rcpp::export]]
std::vector<std::string> cpp_tumopp(const std::vector<std::string>& args) {
    tumopp::Simulation sim(args);
    sim.run();
    return {sim.conf(), sim.tissue().specimens()};
}
