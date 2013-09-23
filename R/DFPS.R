#' Double-Window Flexible Pace Search (DFPS) threshold determination
#'
#' @export
#' @param train_polys training polygons of change areas surrounded by
#' no-change
#' @param chg_img change magnitude image from \code{CVAPS}
#' @param radius radius of no-change surrounding change area polygons
#' @references Chen, J., P. Gong, C. He, R. Pu, and P. Shi. 2003.  
#' Land-use/land-cover change detection using improved change-vector analysis.  
#' Photogrammetric Engineering and Remote Sensing 69:369–380.
#' 
#' Chen, J., X. Chen, X. Cui, and J. Chen. 2011. Change vector analysis in 
#' posterior probability space: a new method for land cover change detection.  
#' IEEE Geoscience and Remote Sensing Letters 8:317–321.
DFPS <- function(train_polys, chg_img, radius=150) {
}
