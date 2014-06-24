#include <RcppArmadillo.h>

using namespace arma;

//' Threshold an image using Huang's fuzzy thresholding method.
//'
//' Implements Huang's fuzzy thresholding method. This function is called by 
//' the \code{\link{threshold}} function. It is not intended to be used 
//' directly.
//'
//' Ported to C++ by from the code in the Auto_threshold imageJ plugin by 
//' Gabriel Landini.
//'
//' See original code at:
//' http://www.mecourse.com/landinig/software/autothreshold/autothreshold.html
//'
//' @param data the input image
//' @return integer threshold value
//' @references Huang, L.-K., and M.-J. J. Wang. 1995. Image thresholding by 
//' minimizing the measures of fuzziness. Pattern recognition 28 (1):41--51.
// [[Rcpp::export]]
int threshold_Huang(arma::ivec& data) {
    // Ported to C++ by Alex Zvoleff from the Auto_threshold imageJ plugin by 
    // Gabriel Landini.
    //
    // See original code at:
    // http://www.mecourse.com/landinig/software/autothreshold/autothreshold.html
    
    // Find first and last non-empty bin
    int first, last;
    for (first = 0; (first < data.n_elem) && (data[first] == 0); first++)
        ; // do nothing
    for (last = data.n_elem - 1; (last > first) && (data[last] == 0); last--)
        ; // do nothing
    if (first == last) return 0;

    //Rcpp::Rcout << "First: " << first << ", Last: " << last << std::endl;
    
    // Calculate the cumulative density and the weighted cumulative density
    fvec S(last + 1);
    fvec W(last + 1);
    S(0) = data(0);
    for (int i = std::max(1, first); i <= last; i++) {
        S(i) = S(i - 1) + data(i);
        W(i) = W(i - 1) + i * data(i);
    }
    
    // Precalculate the summands of the entropy given the absolute difference
    // x - mu (integral)
    double C = last - first;
    fvec Smu(last + 1 - first);
    for (int i = 1; i < Smu.n_elem; i++) {
        double mu = 1 / (1 + abs(i) / C);
        Smu(i) = -mu * log(mu) - (1 - mu) * log(1 - mu);
    }
    
    // Calculate the threshold
    int bestThreshold = 0;
    double bestEntropy = FLT_MAX;
    for (int threshold = first; threshold <= last; threshold++) {
        double entropy = 0;
        int mu = round(W(threshold) / S(threshold));
        for (int i = first; i <= threshold; i++)
            entropy += Smu(abs(i - mu)) * data(i);
        mu = round((W(last) - W(threshold)) / (S(last) - S(threshold)));
        for (int i = threshold + 1; i <= last; i++)
            entropy += Smu(abs(i - mu)) * data(i);

        if (bestEntropy > entropy) {
           bestEntropy = entropy;
           bestThreshold = threshold;
        }
        //Rcpp::Rcout << "debug" << std::endl;
    }
    
    return bestThreshold;
}
