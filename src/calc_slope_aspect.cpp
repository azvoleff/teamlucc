#include <RcppArmadillo.h>
using namespace arma;

//' Calculates slope and aspect using same formula as in the "landsat"
//' package by Sarah Goslee.
//'
//' This function is called by the \code{\link{slopeasp_par}} function. It is 
//' not intended to be used directly.
//'
//' @export
//' @param Rast a matrix of elevations. The slope and aspect will be calculated 
//' for the center pixel in the matrix.
//' @param EWkernel kernel to use for East-West gradient calculation
//' @param EWres East-West resolution
//' @param NSkernel kernel to use for North-South gradient calculation
//' @param NSres North-South resolution
//' @param smoothing positive integer for smoothing. 1 means no smoothing.
//' @return list with two elements. The first element is the calculated slope, 
//' the second element is the calculated aspect.
//' @references
//' Sarah Goslee. Analyzing Remote Sensing Data in {R}: The {landsat} Package.  
//' Journal of Statistical Software, 2011, 43:4, pg 1--25.  
//' http://www.jstatsoft.org/v43/i04/
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
