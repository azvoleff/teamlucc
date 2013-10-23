#include <RcppArmadillo.h>
using namespace arma;

//' Calculates slope and aspect using the same formula as in the "landsat"
//' package by Sarah Goslee.
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector calc_slope_aspect(arma::mat Rast,
        arma::mat EWkernel, double EWres, arma::mat NSkernel,
        double NSres, double smoothing) {
    double EW_mat, NS_mat;
    Rcpp::NumericVector results(4);
    // Calculate the gradient
    EW_mat = accu(Rast % EWkernel) / EWres;
    NS_mat = accu(Rast % NSkernel) / NSres;
    // Calculate the slope (returned as first element of two element vector)
    results[0] = (180.0 / arma::datum::pi) * atan(sqrt(pow(EW_mat, 2) + 
                pow(NS_mat, 2)) / smoothing);
    // Calculate the aspect (returned as second element of two element vector)
    results[1] = 180.0 - (180.0 / arma::datum::pi) * atan(NS_mat / EW_mat) + 
                90.0 * (EW_mat / std::abs(EW_mat));
    if (results[0] == 0)
        results[1] = 0; // set aspect to 0 if slope is 0
    return(results);
}
