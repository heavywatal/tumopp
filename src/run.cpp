// [[Rcpp::plugins(cpp14)]]
#include <Rcpp.h>
#include <tumopp/simulation.hpp>

//' First example
//' @param args string vector
//' @return string
// [[Rcpp::export]]
std::string cpp_tumopp(const std::vector<std::string>& args) {
    tumopp::Simulation sim(args);
    sim.run();
    return sim.tissue().specimens();
}
