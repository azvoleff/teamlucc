#include <RcppArmadillo.h>

using namespace arma;

//' Calculate change direction
//'
//' This code calculate the change direction from two probability images. Not ' 
//' intended to be called directly - see \code{chg_dir}.
//'
//' @export
//' @param t1p time 1 posterior probability matrix (with pixels in rows, bands 
//' in columns)
//' @param t2p time 2 posterior probability matrix (with pixels in rows, bands 
//' in columns)
//' @return vector of change directions
//' @references Chen, J., P. Gong, C. He, R. Pu, and P. Shi. 2003.
//' Land-use/land-cover change detection using improved change-vector analysis.
//' Photogrammetric Engineering and Remote Sensing 69:369-380.
//' 
//' Chen, J., X. Chen, X. Cui, and J. Chen. 2011. Change vector analysis in 
//' posterior probability space: a new method for land cover change detection.  
//' IEEE Geoscience and Remote Sensing Letters 8:317-321.
// [[Rcpp::export]]
arma::ivec calc_chg_dir(arma::mat t1p, arma::mat t2p, int n_classes) {
    ivec chg_dir(t1p.n_rows);
    for (int i = 0; i < t1p.n_rows; i++) {
        mat Eab = (t2p.row(i) - t1p.row(i)) * eye(n_classes, n_classes);
        uword max_location;
        Eab.max(max_location);
        chg_dir(i) = max_location;
    }
    return(chg_dir);
}
