#include <RcppArmadillo.h>
using namespace Rcpp;

bool inC(std::string x, CharacterVector y) {
    bool is_in = false;
    for (int i=0; i < y.size(); ++i) {
        std::string cy = Rcpp::as<std::string>(y(i));
        if (cy == x) {
            is_in = true;
            break;
        }
    }
    return is_in;
}

//' Calculates a glcm texture for use in the glcm.R script
//' @export
// [[Rcpp::export]]
NumericVector calc_texture(NumericMatrix rast, CharacterVector statistics, 
        IntegerVector base_indices, IntegerVector offset_indices,
        int n_grey) {
    arma::mat Rast(rast.begin(), rast.nrow(), rast.ncol(), false); 
    arma::mat G(n_grey, n_grey, arma::fill::zeros);
    arma::mat Pij(n_grey, n_grey, arma::fill::zeros);
    arma::mat imat(n_grey, n_grey, arma::fill::zeros);
    arma::mat jmat(n_grey, n_grey, arma::fill::zeros);
    NumericVector rowsum(n_grey), colsum(n_grey), textures(statistics.size());
    int textures_index=0, Gsum=0;
    double mr=0, mc=0, sig2c=0, sig2r=0;

    for(int i=0; i < offset_indices.size(); i++) {
        G(Rast(base_indices(i)), Rast(offset_indices(i)))++;
    }

    // Calculate Pij from G matrix by dividing each Gij value by number of 
    // total co-ocurrences. Also, make a matrix of i's and a matrix of j's to 
    // be used in the below matrix calculations. These matrices are the same 
    // shape as Pij with the entries equal to the i indices of each cell (for 
    // the imat matrix, which is indexed over the rows) or the j indices of 
    // each cell (for the jmat matrix, which is indexed over the columns).
    Gsum = arma::accu(G);
    for(unsigned i=0; i < G.n_rows; i++) {
        for(unsigned j=0; j < G.n_cols; j++) {
            imat(i, j) = i;
            jmat(i, j) = j;
            Pij(i, j) = G(i, j) / Gsum;
        }
    }

    rowsum = arma::sum(G, 1);
    colsum = arma::sum(G, 0);

    // Calculate mr and mc (forms of col and row means)
    for(unsigned i=0; i < G.n_rows; i++) {
        mr += i * rowsum(i);
        mc += i * colsum(i);
    }

    if (inC("mean", statistics)) {
        // Defined as in Lu and Batistella, 2005, page 252
        textures(textures_index) = mr;
        textures_index++;
    }
    if (inC("variance", statistics)) {
        // Defined as in Haralick, 1973, page 619 (equation 4)
        textures(textures_index) = sum(arma::sum(pow((imat - mr), 2) % Pij, 1));
        textures_index++;
    }
    if (inC("covariance", statistics)) {
        // Defined as in Pratt, 2007, page 540
        textures(textures_index) = sum(arma::sum((imat - mr) % (jmat - mc) % Pij, 1));
        textures_index++;
    }
    if (inC("homogeneity", statistics)) {
        // Defined as in Gonzalez and Woods, 2009, page 832
        textures(textures_index) = sum(arma::sum(Pij / (1 + abs(imat - jmat)), 1));
        textures_index++;
    }
    if (inC("contrast", statistics)) {
        // Defined as in Gonzalez and Woods, 2009, page 832
        textures(textures_index) = sum(arma::sum(pow((imat - jmat), 2) % Pij, 1));
        textures_index++;
    }
    if (inC("dissimilarity", statistics)) {
        //TODO: Find source for dissimilarity
        textures(textures_index) = sum(arma::sum(Pij % abs(imat - jmat), 1));
        textures_index++;
    }
    if (inC("entropy", statistics)) {
        // Defined as in Haralick, 1973, page 619 (equation 9)
        textures(textures_index) = -sum(arma::sum(Pij % log(Pij + .001), 1));
        textures_index++;
    }
    if (inC("second_moment", statistics)) {
        // Defined as in Haralick, 1973, page 619
        textures(textures_index) = sum(arma::sum(pow(Pij, 2), 1));
        textures_index++;
    }
    if (inC("correlation", statistics)) {
        // Defined as in Gonzalez and Woods, 2009, page 832
        for(unsigned i=0; i < G.n_rows; i++) {
            // Calculate sig2r and sig2c (measures of row and column variance)
            sig2r += pow((i - mr), 2) * rowsum(i);
            sig2c += pow((i - mc), 2) * colsum(i);
        }
        textures(textures_index) = sum(arma::sum(((imat - mr) % (jmat - mc) % Pij) /
                    (sqrt(sig2r) * sqrt(sig2c)), 1));
        textures_index++;
    }
    return(textures);
}
