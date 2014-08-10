#include <RcppArmadillo.h>

using namespace arma;

//' Calculate change direction
//'
//' This code calculate the change direction from two probability images. Not
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
arma::ivec calc_chg_dir(arma::mat t1p, arma::mat t2p) {
    ivec chg_dir(t1p.n_rows);
    mat E = eye(t1p.n_cols, t1p.n_cols);
    // the minus below exclues rows where Ea == Eb (persistence)
    mat dEab(E.n_rows * E.n_rows - E.n_rows, t1p.n_cols);
    vec traj_codes(dEab.n_rows);
    int dEab_row = 0;
    // Loop over time 0 (i is t0)
    for (int i = 0; i < E.n_rows; i++) {
        // Loop over time 1 (j is t1)
        for (int j = 0; j < E.n_rows; j++) {
            if (i == j) continue;
            dEab.row(dEab_row) = E.row(j) - E.row(i);
            // Code from=to trajectories by summing t0 and t1 codes after 
            // multiplying t1 codes by the number of classes. Note that 
            // t1p.n_cols is equal to the number of classes.
            traj_codes(dEab_row) = i + j * t1p.n_cols;
            dEab_row++;
        }
    }
    for (int pix_num = 0; pix_num < t1p.n_rows; pix_num++) {
        rowvec dot_dP_dEab(dEab.n_rows);
        rowvec dP = t2p.row(pix_num) - t1p.row(pix_num);
        for (int j = 0; j < dEab.n_rows; j++) {
            dot_dP_dEab(j) = dot(dP, dEab.row(j));
        }
        uword max_location;
        dot_dP_dEab.max(max_location);
        chg_dir(pix_num) = traj_codes(max_location);
    }
    return(chg_dir);
}
