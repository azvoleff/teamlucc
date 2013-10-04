#' Fractional mixing model implemented with support vector machine (SVM)
#'
#' This code calculate the change magnitude image for the Change Vector 
#' Analysis in Posterior Probability Space (CVAPS) method of Chen et al. 2011.  
#' Use the change magnitude image in conjunction with the change direction 
#' image from \code{chg_dir}, and \code{DFPS} to use the Double Window Flexible 
#' Pace Search method (Chen et al. 2003) to determine the threshold to use to 
#' map areas of change and no-change.
#'
#' @export
#' @param x a \code{Raster*} image with the predictor layer(s) for the 
#' classification
#' @param train_data a data table with a column labeled 'y' with the observed 
#' classes, and one or more columns with the values of predictor(s) at each 
#' location.
#' @return a \code{RasterBrick} with the proportion of each pixel within each 
#' of the classes listed in the training data
#' @param filename (optional) filename for output image
#' @references Brown, M., H. G.  Lewis, and S.  R.  Gunn.  2000.  Linear 
#' spectral mixture models and support vector machines for remote sensing.  
#' Geoscience and Remote Sensing, IEEE Transactions on 38:2346–2360.
mixing_model <- function(x, train_data, filename=NULL) {
    #TODO: Code this function
}
