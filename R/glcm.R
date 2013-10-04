#' Image texture measures from grey-level co-occurrence matrices (GLCM)
#'
#' @export
#' @param x a /code{RasterLayer}
#' @param n_grey number of grey levels to use in texture calculation
#' @param window the window size to consider for texture calculation as a two 
#' element integer vector
#' @param shift a 2 element integer vector giving the shift (Q in Gonzalez and 
#' Woods, 2008)
#' @param statistic A list of GLCM texture measures to calculate.  Can include 
#' any (one or more) of the following: 'mean', 'variance', 'homogeneity', 
#' 'contrast', 'dissimilarity', 'entropy', 'second_moment', 'correlation'.
#' @return A /code{RasterLayer} with the calculated requested GLCM texture 
#' measures.
#' @references Lu, D., and M. Batistella. 2005. Exploring TM image texture and 
#' its relationships with biomass estimation in Rondônia, Brazilian Amazon.  
#' Acta Amazonica 35:249–257.
#'
#' Gonzalez, R. C. 2008. Digital image processing. 3rd ed. Prentice Hall, Upper 
#' Saddle River, N.J, page 830-836.
#'
#' Haralick, R. M., K. Shanmugam, and I. Dinstein. 1973. Textural features for 
#' image classification. IEEE Transactions on Systems, Man and Cybernetics 
#' SMC-3:610 –621.
#' @examples
library(teamr)
L5TSR_1986 <- stack(system.file('extdata/L5TSR_1986.dat', package='teamr'))
x <- raster(L5TSR_1986, layer=1)
n_grey <- 3
window <- c(3, 3)
shift <- c(1, 1)
statistic <- 'mean'
x <- raster(matrix(c(0, 0, 0, 0, 0,
              0, 0, 0, 2, 0,
              0, 0, 2, 2, 0,
              1, 1, 2, 3, 0,
              1, 1, 2, 3, 0), nrow=5))


glcm <- function(x, n_grey=32, window=c(3, 3), shift=c(1, 1),
                 statistic='mean') {
    # Define shift[2] to be in the North (up) direction:
    shift[2] <- -shift[2]

    if (length(window) != 2) {
        stop('window must be integer vector of length 2')
    }
    if (length(shift) != 2) {
        stop('shift must be integer vector of length 2')
    }
    if ((window[1] < 3) || (window[2] < 3)) {
        stop('both elements of window must be  >= 3')
    }
    if ((window[1] %% 2 == 0) || (window[2] %% 2 == 0)) {
        stop('both elements of window must be odd')
    }

    # Resample the image to the required number of grey levels
    x_grey <- cut(x, breaks=seq(cellStats(x, 'min'), cellStats(x, 'max'), 
                                length.out=n_grey + 1), include.lowest=TRUE)

    rast <- matrix(getValues(x_grey), nrow=5)[c(1:3), c(2:4)]

    # Setup coordinates for offset image
    col_coords <- matrix(rep(shift[1] + seq(1, window[2]), window[1]), 
                       nrow=window[1], byrow=TRUE)
    row_coords <- matrix(rep(shift[2] + seq(1, window[1]), window[2]), 
                       nrow=window[1])
    # Set coordinates that fall outside the window to NA
    col_coords[col_coords > window[2] | col_coords < 1] <- NA
    row_coords[row_coords > window[1] | row_coords < 1] <- NA
    col_coords[is.na(row_coords)] <- NA
    row_coords[is.na(col_coords)] <- NA

    calc_texture <- function(rast, statistics, col_coords, row_coords, ...) {
        G <- matrix(0, n_grey, n_grey)
        co_occur <- cbind(as.vector(rast), rast[cbind(as.vector(row_coords), 
                                                      as.vector(col_coords))])
        for (n in 1:nrow(co_occur)) {
            G[co_occur[n, 1], co_occur[n, 2]] <- G[co_occur[n, 1], co_occur[n, 2]] + 1
        }
        Pij = G / sum(G)

        # Calculate rowSums and colSums of Pij
        rowsum = rowSums(Pij)
        colsum = colSums(Pij)
        # Calcuate mr and mc (forms of col and row means) and sig2r and sig2c 
        # (measures of row and column variance)
        mr = sum(c(1:nrow(Pij)) * rowsum)
        mc = sum(c(1:ncol(Pij)) * colsum)
        sig2r = sum((c(1:nrow(Pij)) - mr)^2 * rowsum)
        sig2c = sum((c(1:nrow(Pij)) - mc)^2 * colsum)

        if ('mean' %in% statistics) {
        }
        if ('variance' %in% statistics) {
        }
        if ('homogeneity' %in% statistics) {
        }
        if ('contrast' %in% statistics) {
        }
        if ('dissimilarity' %in% statistics) {
        }
        if ('entropy' %in% statistics) {
        }
        if ('second_moment' %in% statistics) {
        }
        if ('correlation' %in% statistics) {
        }
    }
    textures <- rasterEngine(x_grey, fun=calc_texture,
                             args=list(statistics, col_coords, row_coords), window_dims=window)
    
    return(textures)
}
a  <- matrix(seq(1, 25), nrow=5)
a[c(3:5),c(3:5)]
# Encode matrix from Fig. 3 in Haralick (1973) for testing:
#a <- matrix(c(0, 0, 0, 2, 0, 0, 2, 2, 1, 1, 2, 3, 1, 1, 2, 3), nrow=4)
a <- matrix(c(0, 0, 0, 0, 0,
              0, 0, 0, 2, 0,
              0, 0, 2, 2, 0,
              1, 1, 2, 3, 0,
              1, 1, 2, 3, 0), nrow=5)

# n_grey is the number of grey levels in the image
n_grey <- 4
# d is distance value (in pixels)
d <- 2
# w is number of pixels to consider in each direction (horizontal and 
# vertical), and including the center pixel, but not taking into account d (for 
# this see 's' below).
w <- 3
# s is the total number of pixels needed in each direction (vertical upwards 
# and horizontal to the right) for all comparisons to be valid and to not fall 
# outside of the original image matrix.
s <- (w + d)

# o_x is origin x value (center x coordinate)
o_x <- 3
# o_y is origin y value (center y coordinate)
o_y <- 3

# First make a matrix with the coordinates of all the pixels in the window, 
# relative to the coordinates of the central pixel, considered to be at (0, 0).
w_x <- matrix(rep(seq(-w, w), s), nrow=s, byrow=TRUE)
w_y <- matrix(rep(seq(-w, w), s), nrow=s, byrow=FALSE)

# Now make a matrix with the coordinates for a 0 degree neighbor, the neighbor 
# immediately to the left of each pixel.
nb_x_0 <- w_x + 1
nb_y_0 <- w_y

w_vals <- a[cbind(as.vector(w_y + o_y), as.vector(w_x + o_x))]
nb_0 <- a[cbind(as.vector(nb_y_0 + o_y), as.vector(nb_x_0 + o_x))]
nb_45 <- a[cbind(as.vector(nb_y_45 + o_y), as.vector(nb_x_45 + o_x))]
nb_90 <- a[cbind(as.vector(nb_y_90 + o_y), as.vector(nb_x_90 + o_x))]
nb_135 <- a[cbind(as.vector(nb_y_135 + o_y), as.vector(nb_x_135 + o_x))]

occ[w_vals, nb_0] <- occ[w_vals, nb_0] + 1
occ[w_vals, nb_45] <- occ[w_vals, nb_45] + 1
occ[w_vals, nb_90] <- occ[w_vals, nb_90] + 1
occ[w_vals, nb_135] <- occ[w_vals, nb_135] + 1

# Make symmetric to account for forward and backward relationships:
occ[upper.tri(occ)] <- occ[upper.tri(occ)] + t(occ)[upper.tri(occ)]
occ[lower.tri(occ)] <- t(occ)[lower.tri(occ)]
occ[diag(nrow(occ))] <- 2 * diag(occ)
