// [[Rcpp::plugins(cpp14)]]
#include <Rcpp.h>
#include <tumopp/simulation.hpp>

//' First example
//' @param x numeric vector
//' @return numeric vector
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector tumopp_test(Rcpp::NumericVector x) {
    tumopp::Simulation sim(0, {});
    sim.run();
    return x;
}
