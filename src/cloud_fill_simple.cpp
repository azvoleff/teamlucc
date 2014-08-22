#include <RcppArmadillo.h>

using namespace arma;

//' Cloud fill using a simple linear model approach
//'
//' This algorithm fills clouds using a simple approach in which the value of 
//' each clouded pixel is calculated using a linear model. The script
//' develops a separate linear model (with slope and intercept) for each band 
//' and each cloud. For each cloud, and each image band, the script finds all 
//' pixels clear in both the cloudy and fill images, and calculates a 
//' regression model in which pixel values in the fill image are the 
//' independent variable, and pixel values in the clouded image are the 
//' dependent variable. The script then uses this model to predict pixel values 
//' for each band in each cloud in the clouded image.
//'
//' This function is called by the \code{\link{cloud_remove}} function. It is
//' not intended to be used directly.
//'
//' @param cloudy the cloudy image as a matrix, with pixels in columns (in 
//' column-major order) and with number of columns equal to number of bands
//' @param clear the clear image as a matrix, with pixels in columns (in 
//' column-major order) and with number of columns equal to number of bands
//' @param cloud_mask the cloud mask image as a vector (in column-major order), 
//' with clouds coded with unique integer codes starting at 1, and with areas 
//' that are clear in both images  coded as 0. Areas that are missing in the 
//' clear image, should be coded as -1.
//' @param dims the dimensions of the cloudy image as a length 3 vector: (rows, 
//' columns, bands)
//' @param num_class set the estimated number of classes in image
//' @param cloud_nbh the range of cloud neighborhood (in pixels)
//' @param DN_min the minimum valid DN value
//' @param DN_max the maximum valid DN value
//' @param verbose whether to print detailed status messages
//' @return array with cloud filled image with dims: cols, rows, bands
//' parameter, containing the selected textures measures
// [[Rcpp::export]]
arma::mat cloud_fill_simple(arma::mat cloudy, arma::mat& clear,
        arma::ivec& cloud_mask, arma::ivec dims, int num_class,
        int cloud_nbh, int DN_min, int DN_max, bool verbose=false) {

    if (verbose) Rcpp::Rcout << "dims: " << dims(0) << ", " << dims(1) << ", " 
        << dims(2) << std::endl;
    if (verbose) Rcpp::Rcout << "cloudy dims: " << cloudy.n_rows << ", " << 
        cloudy.n_cols << std::endl;

    // Make a list of the cloud codes in this file - anything less than 1 is 
    // not a cloud code (0 is no clear, and -1 means no data in the clear 
    // image)
    ivec cloud_codes = unique(cloud_mask);
    cloud_codes = cloud_codes(find(cloud_codes >= 1));

    // Allow treating the multiband images as cubes
    cube cloudy_cube(cloudy.begin(), dims(0), dims(1), dims(2), false);
    cube clear_cube(clear.begin(), dims(0), dims(1), dims(2), false);
    // Allow also treating cloud_mask as 2d matrix (row, cols)
    imat cloud_mask_mat(cloud_mask.begin(), dims(0), dims(1), false);

    if (verbose) Rcpp::Rcout << cloud_codes.n_elem  << " cloud(s) to fill" << std::endl;

    for (unsigned n=0; n < cloud_codes.n_elem; n++) {
        int cloud_code = cloud_codes(n);
        if (verbose) Rcpp::Rcout << "Filling cloud " << cloud_code;

        // These indices refer to the position of cloud pixels within the 
        // overall cloud_mask block
        uvec cloud_vec_i = find(cloud_mask == cloud_code);
        uvec cloud_col_i = floor(cloud_vec_i / dims(0));
        uvec cloud_row_i = cloud_vec_i - cloud_col_i * dims(0);

        // Now add in the cloud neighborhood
        int left_col = min(cloud_col_i) - cloud_nbh;
        if (left_col < 0) left_col = 0;

        int right_col = max(cloud_col_i) + cloud_nbh;
        if (right_col > (dims(1) - 1)) right_col = dims(1) - 1;

        int up_row = min(cloud_row_i) - cloud_nbh;
        if (up_row < 0) up_row = 0;

        int down_row = max(cloud_row_i) + cloud_nbh;
        if (down_row > (dims(0) - 1)) down_row = dims(0) - 1;

        int num_sub_cols = (right_col - left_col) + 1;
        int num_sub_rows = (down_row - up_row) + 1;
        // Extract the cloud neighborhood from the cubes, and setup 
        // column-major matrices that will be used for the remaining 
        // calculations.
        cube sub_cloudy_cube = cloudy_cube.tube(up_row, left_col, down_row, right_col);
        cube sub_clear_cube = clear_cube.tube(up_row, left_col, down_row, right_col);
        mat sub_cloudy(num_sub_rows*num_sub_cols, dims(2));
        mat sub_clear(num_sub_rows*num_sub_cols, dims(2));
        for (unsigned elnum=0; elnum < sub_clear_cube.n_elem; elnum++) {
            sub_cloudy(elnum) = sub_cloudy_cube(elnum);
            sub_clear(elnum) = sub_clear_cube(elnum);
        }

        imat sub_cloud_mask = cloud_mask_mat.submat(up_row, left_col, down_row, right_col);

        // These indices refer to the position of pixels of this cloud
        // within the cloud neighborhood of this cloud (a subset of cloud_mask 
        // block)
        uvec sub_cloud_vec_i = find(sub_cloud_mask == cloud_code);
        uvec sub_cloud_col_i = floor(sub_cloud_vec_i / sub_cloud_mask.n_rows);
        uvec sub_cloud_row_i = sub_cloud_vec_i - sub_cloud_col_i * sub_cloud_mask.n_rows;

        if (verbose) Rcpp::Rcout << " (" << sub_cloud_vec_i.n_elem <<  " pixels)" << std::endl;

        // These indices refer to the position of clear pixels within the cloud 
        // neighborhood of this cloud (a subset of clear block)
        uvec sub_clear_vec_i = find(sub_cloud_mask == 0);
        if (sub_clear_vec_i.n_elem == 0) {
            if (verbose) Rcpp::Rcout << "No clear neighbors in cloudy image. Skipping fill." << std::endl;
            continue;
        }
        uvec sub_clear_col_i = floor(sub_clear_vec_i / sub_cloud_mask.n_rows);
        uvec sub_clear_row_i = sub_clear_vec_i - sub_clear_col_i * sub_cloud_mask.n_rows;

        mat sub_clear_clear = sub_clear.rows(sub_clear_vec_i);
        mat sub_cloudy_clear = sub_cloudy.rows(sub_clear_vec_i);

        mat sub_clear_cloudy = sub_clear.rows(sub_cloud_vec_i);

        // lm code is based on code from Dirk Eddelbuettel at: 
        // http://bit.ly/1oTa60F
        for(int iband=0; iband < dims(2); iband++) {
            mat X_param = join_rows(sub_clear_clear.col(iband), ones(sub_clear_clear.n_rows));
            colvec coef;
            try {
                coef = solve(X_param, sub_cloudy_clear.col(iband)); // fit model y ~ X + 1
            } catch(std::exception &ex) {	
                // cannot solve (singular), so assume slope 1, intercept 0
                if (verbose) Rcpp::Rcout << "solve() failed - assuming slope 1, intercept 0." << std::endl;
                coef << 1 << endr << 0 << endr;
            } catch(...) { 
                ::Rf_error("c++ exception (unknown reason)"); 
            }

            mat X_pred = join_rows(sub_clear_cloudy.col(iband), ones(sub_clear_cloudy.n_rows));
            colvec preds = X_pred * coef; // make predictions

            for (unsigned ic=0; ic < sub_cloud_vec_i.n_elem; ic++) {
                // Calculate row and column location of target pixel
                cloudy_cube(up_row + sub_cloud_row_i(ic), left_col + sub_cloud_col_i(ic), iband) = preds(ic);
            }
            // cloudy_cube.tube(up_row, left_col, down_row, right_col).elem(sub_cloud_vec_i) = preds;
        }

        if (verbose) Rcpp::Rcout << std::endl;
    }
    return(cloudy);
}
