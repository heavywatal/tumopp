// [[Rcpp::plugins(cpp14)]]
#include <Rcpp.h>
#include <tumopp/simulation.hpp>

//' First example
//' @param args numeric vector
//' @return numeric vector
//' @export
// [[Rcpp::export]]
std::string run_tumopp(Rcpp::CharacterVector args=Rcpp::CharacterVector::create()) {
    auto vs_args = Rcpp::as<std::vector<std::string>>(args);
    tumopp::Simulation sim(vs_args);
    sim.run();
    return sim.tissue().specimens();
}
