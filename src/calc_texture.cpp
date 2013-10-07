#include <Rcpp.h>
using namespace Rcpp;

// Adapted from http://adv-r.had.co.nz/Rcpp.html
NumericVector colSumsC(NumericMatrix x) {
    int ncol = x.ncol(), nrow = x.nrow();
    NumericVector out(ncol);
    for (int i = 0; i < ncol; i++) {
        double total = 0;
        for (int j = 0; j < nrow; j++) {
            total += x(i, j);
        }
        out[i] = total;
    }
    return out;
}

// From http://adv-r.had.co.nz/Rcpp.html
NumericVector rowSumsC(NumericMatrix x) {
    int nrow = x.nrow(), ncol = x.ncol();
    NumericVector out(nrow);
    for (int i = 0; i < nrow; i++) {
        double total = 0;
        for (int j = 0; j < ncol; j++) {
            total += x(i, j);
        }
        out[i] = total;
    }
    return out;
}


// [[Rcpp::export]]
NumericVector calc_textureC(NumericMatrix rast, CharacterVector statistics, 
        IntegerVector base_indices, IntegerVector offset_indices,
        Integer n_grey) {
    NumericMatrix G(n_grey, n_grey), imat=(n_grey, n_grey),
    jmat=(n_grey, n_grey);
    int base_value = 0, offset_value = 0;
    NumericVector rowsum(n_grey), colsum(n_grey), textures(statistics.size());
    int textures_index = 0;
    double mr, mc, sig2c, sig2r;

    for(int i=0; i < offset_indices.size(); i++) {
        base_value = rast[,,1][base_indices];
        offset_value = rast[,,1][offset_indices];
        G(base_value, offset_value)++;
    }
    Pij = G / sum(G);

    // Calculate rowSums and colSums of Pij
    rowsum = rowSumsC(G);
    colsum = colSumsC(G);

    // Calculate mr and mc (forms of col and row means) and sig2r and sig2c 
    // (measures of row and column variance)
    mr = sum(c(1:nrow(Pij)) * rowsum);
    mc = sum(c(1:ncol(Pij)) * colsum);
    sig2r = sum((c(1:nrow(Pij)) - mr)^2 * rowsum);
    sig2c = sum((c(1:ncol(Pij)) - mc)^2 * colsum);

    // Make a matrix of i's and a matrix of j's to be used in the below 
    // matrix calculations. These matrices are the same shape as Pij with 
    // the entries equal to the i indices of each cell (for the imat matrix, 
    // which is indexed over the rows) or the j indices of each cell (for 
    // the jmat matrix, which is indexed over the columns).
    imat = matrix(rep(1:nrow(Pij), ncol(Pij)), nrow=nrow(Pij));
    jmat = matrix(rep(1:ncol(Pij), nrow(Pij)), ncol=ncol(Pij), byrow=TRUE);

    if ("mean" %in% statistics) {
        // Defined as in Lu and Batistella, 2005, page 252
        textures(textures_index) = mr;
        textures_index ++;
    }
    if ("variance" %in% statistics) {
        // Defined as in Haralick, 1973, page 619 (equation 4)
        textures(textures_index) = sum(rowSums((imat - mr)^2 * Pij));
        textures_index ++;
    }
    if ("covariance" %in% statistics) {
        // Defined as in Pratt, 2007, page 540
        textures(textures_index) = sum(rowSums((imat - mr) *
                    (jmat - mc) * Pij));
        textures_index ++;
    }
    if ("homogeneity" %in% statistics) {
        // Defined as in Gonzalez and Woods, 2009, page 832
        textures(textures_index) = sum(rowSums(Pij / (1 + abs(imat - jmat))));
        textures_index ++;
    }
    if ("contrast" %in% statistics) {
        // Defined as in Gonzalez and Woods, 2009, page 832
        textures(textures_index) = sum(rowSums((imat - jmat)^2 * Pij));
        textures_index ++;
    }
    if ("dissimilarity" %in% statistics) {
        //TODO: Find source for dissimilarity
        textures(textures_index) = sum(rowSums(Pij * abs(imat - jmat)));
        textures_index ++;
    }
    if ("entropy" %in% statistics) {
        // Defined as in Haralick, 1973, page 619 (equation 9)
        textures(textures_index) = -sum(rowSums(Pij * log(Pij + .001)));
        textures_index ++;
    }
    if ("second_moment" %in% statistics) {
        // Defined as in Haralick, 1973, page 619
        textures(textures_index) = sum(rowSums(Pij^2));
        textures_index ++;
    }
    if ("correlation" %in% statistics) {
        // Defined as in Gonzalez and Woods, 2009, page 832
        textures(textures_index) = sum(rowSums(((imat - mr) * (jmat - mc) * Pij) /
                    (sqrt(sig2r) * sqrt(sig2c))));
        textures_index ++;
    }
}
