// [[Rcpp::plugins(cpp14)]]
#include <Rcpp.h>
#include <tumopp/simulation.hpp>

//' Run C++ simulation
//' @return conf and population as strings
//' @rdname tumopp
// [[Rcpp::export]]
std::vector<std::string> cpp_tumopp(const std::vector<std::string>& args, const size_t npair=0) {
    tumopp::Simulation sim(args);
    sim.run();
    return {
        sim.conf(),
        sim.tissue().specimens(),
        sim.tissue().drivers(),
        sim.tissue().pairwise_distance(npair)};
}
