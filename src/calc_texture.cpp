#include <RcppArmadillo.h>
using namespace arma;

// Define a pointer to a texture function that will be used to map selected 
// co-occurrence statistics to the below texture calculation functions.
typedef double (*pfunc)(mat, mat, mat, double);

double text_mean(mat pij, mat imat, mat jmat, double mu) {
    // Defined as in Lu and Batistella, 2005, page 252
    return(mu);
}

double text_variance(mat pij, mat imat, mat jmat, double mu) {
    // Defined as in Haralick, 1973, page 619 (equation 4)
    return(accu(square(imat - mu) % pij));
}

double text_covariance(mat pij, mat imat, mat jmat, double mu) {
    // Defined as in Pratt, 2007, page 540
    return(accu((imat - mean(pij, 1)) % (jmat - mean(pij, 0)) % pij));
}

double text_homogeneity(mat pij, mat imat, mat jmat, double mu) {
    // Defined as in Haralick, 1973, page 619 (equation 5)
    return(accu(pij / (1 + square(imat - jmat))));
}

double text_contrast(mat pij, mat imat, mat jmat, double mu) {
    // Defined as in Gonzalez and Woods, 2009, page 832
    // Note that pij.n_rows is equivalent to n_grey
    // return(accu(square(linspace<vec>(0, pij.n_rows-1, pij.n_rows)) * (imat - 
    //                 jmat % pij)));
    return(accu(square(linspace<vec>(0, pij.n_rows-1, pij.n_rows)) * (imat - jmat % pij)));
}

double text_dissimilarity(mat pij, mat imat, mat jmat, double mu) {
    //TODO: Find source for dissimilarity
    return(accu(pij % abs(imat - jmat)));
}

double text_entropy(mat pij, mat imat, mat jmat, double mu) {
    // Defined as in Haralick, 1973, page 619 (equation 9)
    mat pij_log(pij);
    pij_log = mat(pij);
    pij_log(find(pij_log)) = log(pij_log(find(pij_log)));
    return(-accu(pij % pij_log));
}

double text_second_moment(mat pij, mat imat, mat jmat, double mu) {
    // Defined as in Haralick, 1973, page 619
    return(accu(square(pij)));
}

double text_correlation(mat pij, mat imat, mat jmat, double mu) {
    // Defined as in Gonzalez and Woods, 2009, page 832
    // Calculate mr and mc (forms of col and row means), see Gonzalez and 
    // Woods, 2009, page 832
    double sigc, sigr, mr, mc;
    mr = sum(linspace<vec>(1, pij.n_cols, pij.n_cols) % sum(pij, 1));
    mc = sum(trans(linspace<vec>(1, pij.n_rows, pij.n_rows)) % sum(pij, 0));
    // Calculate sigr and sigc (measures of row and column std deviation)
    sigr = sqrt(sum(square(linspace<vec>(1, pij.n_cols, pij.n_cols) - mr) % sum(pij, 1)));
    sigc = sqrt(sum(square(trans(linspace<vec>(1, pij.n_rows, pij.n_rows)) - mc) % sum(pij, 0)));
    return((accu(imat % jmat % pij) - mr * mc) / (sigr * sigc));
}

//' Calculates a glcm texture for use in the glcm.R script
//'
//' This function is called by the \code{\link{glcm}} function. It is 
//' not intended to be used directly.
//'
//' @export
//' @param rast a matrix containing the pixels to be used in the texture 
//' calculation
//' @param statistics a list of strings naming the texture statistics to 
//' calculate
//' @param n_grey number of grey levels to use in texture calculation
//' @param window_dims 2 element list with row and column dimensions of the
//' texture window
//' @param shift a length 2 vector with the number of cells to shift when
//' computing co-ocurrency matrices
//' @return a list of length equal to the length of the \code{statistics} input 
//' parameter, containing the selected textures measures
//' @references
//' Sarah Goslee. Analyzing Remote Sensing Data in {R}: The {landsat} Package.  
//' Journal of Statistical Software, 2011, 43:4, pg 1--25.  
//' http://www.jstatsoft.org/v43/i04/
// [[Rcpp::export]]
arma::cube calc_texture_full_image(arma::mat rast,
        Rcpp::CharacterVector statistics, int n_grey,
        arma::vec window_dims, arma::vec shift) {
    mat imat(n_grey, n_grey);
    mat jmat(n_grey, n_grey);
    mat base_window(window_dims(0), window_dims(1));
    mat offset_window(window_dims(0), window_dims(1));
    mat pij(n_grey, n_grey);
    vec base_ul(2), offset_ul(2), center_coord(2);
    double mu;
    // textures cube will hold the calculated texture statistics
    cube textures(rast.n_rows, rast.n_cols, statistics.size(), fill::zeros);

    std::map<std::string, double (*)(mat, mat, mat, double)> stat_func_map;
    stat_func_map["mean"] = text_mean;
    stat_func_map["variance"] = text_variance;
    stat_func_map["covariance"] = text_covariance;
    stat_func_map["homogeneity"] = text_homogeneity;
    stat_func_map["contrast"] = text_contrast;
    stat_func_map["dissimilarity"] = text_dissimilarity;
    stat_func_map["entropy"] = text_entropy;
    stat_func_map["second_moment"] = text_second_moment;
    stat_func_map["correlation"] = text_correlation;

    // Calculate the base upper left (ul) coords and offset upper left coords 
    // as row, column with zero based indices.
    base_ul = vec("0 0");
    if (shift[0] < 0) {
        base_ul[0] = base_ul[0] + abs(shift[0]);
    }
    if (shift[1] < 0) {
        base_ul[1] = base_ul[1] + abs(shift[1]);
    }
    offset_ul = base_ul + shift;
    center_coord = base_ul + floor(window_dims / 2);
     
    // Make a matrix of i's and a matrix of j's to be used in the below matrix 
    // calculations. These matrices are the same shape as pij with the entries 
    // equal to the i indices of each cell (for the imat matrix, which is 
    // indexed over the rows) or the j indices of each cell (for the jmat 
    // matrix, which is indexed over the columns). Note that linspace<mat> 
    // makes a column vector.
    imat = repmat(linspace<vec>(1, pij.n_rows, pij.n_rows), 1, pij.n_cols);
    jmat = trans(imat);

    for(unsigned row=0; row < (rast.n_rows - abs(shift(0)) - ceil(window_dims(0)/2)); row++) {
        // if (row %250 == 0 ) {
        //     Rcpp::Rcout << "Row: " << row << std::endl;
        // }
        for(unsigned col=0; col < (rast.n_cols - abs(shift(1)) - ceil(window_dims(1)/2)); col++) {
            base_window = rast.submat(row + base_ul(0),
                                      col + base_ul(1),
                                      row + base_ul(0) + window_dims(0) - 1,
                                      col + base_ul(1) + window_dims(1) - 1);
            offset_window = rast.submat(row + offset_ul(0),
                                        col + offset_ul(1),
                                        row + offset_ul(0) + window_dims(0) - 1,
                                        col + offset_ul(1) + window_dims(1) - 1);
            pij.fill(0);
            for(unsigned i=0; i < base_window.n_elem; i++) {
                // Subtract one from the below indices to correct for row and col 
                // indices starting at 0 in C++ versus 1 in R.
                pij(base_window(i) - 1, offset_window(i) - 1)++;
            }
            pij = pij / base_window.n_elem;

            mu = accu(base_window) / base_window.n_elem;

            // Loop over the selected statistics, using the stat_func_map map 
            // to map each selected statistic to the appropriate texture 
            // function.
            for(signed i=0; i < statistics.size(); i++) {
                pfunc f = stat_func_map[Rcpp::as<std::string>(statistics(i))];
                textures(row + center_coord(0),
                         col + center_coord(1), i) = (*f)(pij, imat, jmat, mu);
            }

        }
    }
    return(textures);
}
