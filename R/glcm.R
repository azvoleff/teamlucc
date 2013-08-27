#' Calculate texture measures from grey-level co-occurrence matrices (GLCM)
#'
#' @export
#' @param img A /code{RasterLayer}.
#' @param n_grey Number of grey levels to use.
#' @param d Distance value (in pixels).
#' @param w The number of pixels to consider in each direction (horizontal and 
#' vertical), and including the center pixel, but not taking into account 
#' /code{d}).
#' @param statistic A list of GLCM texture measures to calculate.
#' @return A /code{RasterLayer} with the calculated requested GLCM texture 
#' measures.
#' @examples
#' # TODO: Need an example
glcm <- function(img, n_grey, d, w, statistic) {
    # TODO: Finish this
    return(1)
}
# a  <- matrix(seq(1, 25), nrow=5)
# a[c(3:5),c(3:5)]
# # Encode matrix from Fig. 3 in Haralick (1973) for testing:
# #a <- matrix(c(0, 0, 0, 2, 0, 0, 2, 2, 1, 1, 2, 3, 1, 1, 2, 3), nrow=4)
# a <- matrix(c(0, 0, 0, 0, 0,
#               0, 0, 0, 2, 0,
#               0, 0, 2, 2, 0,
#               1, 1, 2, 3, 0,
#               1, 1, 2, 3, 0), nrow=5)
# 
# # n_gray is the number of grey levels in the image
# n_gray <- 4
# # d is distance value (in pixels)
# d <- 2
# # w is number of pixels to consider in each direction (horizontal and 
# # vertical), and including the center pixel, but not taking into account d (for 
# # this see 's' below).
# w <- 3
# # s is the total number of pixels needed in each direction (vertical upwards 
# # and horizontal to the right) for all comparisons to be valid and to not fall 
# # outside of the original image matrix.
# s <- (w + d)
# 
# # o_x is origin x value (center x coordinate)
# o_x <- 
# # o_y is origin y value (center y coordinate)
# o_y <- 
# 
# # Check that the data matrix is large enough to support the comparisons:
# if (
# 
# # First make a matrix with the coordinates of all the pixels in the window, 
# # relative to the coordinates of the central pixel, considered to be at (0, 0).
# w_x <- matrix(rep(seq(-w, w), s), nrow=s, byrow=TRUE)
# w_y <- matrix(rep(seq(-w, w), s), nrow=s, byrow=FALSE)
# 
# # Now make a matrix with the coordinates for a 0 degree neighbor, the neighbor 
# # immediately to the left of each pixel.
# nb_x_0 <- w_x + 1
# nb_y_0 <- w_y
# 
# # To rotate 45 degrees, need to subtract 1 from to the y coordinates only (as 
# # the y coordinates decrease (lower numbered rows) going in the vertical 
# # direction that would otherwise be positive with Cartesian coordinates.
# nb_x_45 <- w_x + d
# nb_y_45 <- w_y - d
# nb_x_90 <- w_x
# nb_y_90 <- w_y - d
# nb_x_135 <- w_x - d
# nb_y_135 <- w_y - d
# 
# #w_vals <- matrix(a[cbind(as.vector(w_y + o_y), as.vector(w_x + o_x))], nrow=s)
# #nb_0 <- matrix(a[cbind(as.vector(nb_y_0 + o_y), as.vector(nb_x_0 + o_x))], nrow=s)
# #nb_45 <- matrix(a[cbind(as.vector(nb_y_45 + o_y), as.vector(nb_x_45 + o_x))], nrow=s)
# #nb_90 <- matrix(a[cbind(as.vector(nb_y_90 + o_y), as.vector(nb_x_90 + o_x))], nrow=s)
# #nb_135 <- matrix(a[cbind(as.vector(nb_y_135 + o_y), as.vector(nb_x_135 + o_x))], nrow=s)
# 
# w_vals <- a[cbind(as.vector(w_y + o_y), as.vector(w_x + o_x))]
# nb_0 <- a[cbind(as.vector(nb_y_0 + o_y), as.vector(nb_x_0 + o_x))]
# nb_45 <- a[cbind(as.vector(nb_y_45 + o_y), as.vector(nb_x_45 + o_x))]
# nb_90 <- a[cbind(as.vector(nb_y_90 + o_y), as.vector(nb_x_90 + o_x))]
# nb_135 <- a[cbind(as.vector(nb_y_135 + o_y), as.vector(nb_x_135 + o_x))]
# 
# occ <- matrix(0, n_gray, n_gray)
# occ[w_vals, nb_0] <- occ[w_vals, nb_0] + 1
# occ[w_vals, nb_45] <- occ[w_vals, nb_45] + 1
# occ[w_vals, nb_90] <- occ[w_vals, nb_90] + 1
# occ[w_vals, nb_135] <- occ[w_vals, nb_135] + 1
# 
# # Make symmetric to account for forward and backward relationships:
# occ[upper.tri(occ)] <- occ[upper.tri(occ)] + t(occ)[upper.tri(occ)]
# occ[lower.tri(occ)] <- t(occ)[lower.tri(occ)]
# occ[diag(nrow(occ))] <- 2 * diag(occ)
