#' Draw a random sample from a grid laid out on a RasterLayer or matrix
#'
#' This function is used to subsample a \code{RasterLayer} or \code{matrix} by 
#' dividing the dataset into a grid of \code{horizcells} x \code{vertcells}, 
#' and by then drawing a sample of size \code{nsamp} from within each grid 
#' cell.
#'
#' @export
#' @import rgeos
#' @param x a matrix or RasterLayer to draw sample from
#' @param horizcells how many cells to break the raster in horizontally (over 
#' the columns)
#' @param vertcells how many cells to break the raster in vertically (over 
#' the rows)
#' @param nsamp how many samples to draw from each grid cell
#' @return a vector of sampled values of length \code{horizcells} * 
#' \code{vertcells} * \code{nsamp}
#' @examples
#' # Make a 100x100 matrix
#' x <- matrix(1:10000, nrow=100)
#' # Subsample the matrix by breaking it up into a 10 x 10 grid, and take 10
#' # random samples from each grid cell without replacement (these are the
#' # default parameters).
#' y <- gridsample(x)
gridsample <- function(x, horizcells=10, vertcells=10, nsamp=10) {
    horizstart <- round(seq(1, ncol(x), ncol(x) / horizcells))
    horizend <- c((horizstart -1)[2:length(horizstart)], ncol(x))
    vertstart <- round(seq(1, nrow(x), nrow(x) / vertcells))
    vertend <- c((vertstart -1)[2:length(vertstart)], nrow(x))
    retvals <- vector('numeric', horizcells * vertcells * nsamp)
    retval_index <- 1
    for (horizcellnum in 1:length(horizstart)) {
        for (vertcellnum in 1:length(vertstart)) {
            cell_ncols <- horizend[horizcellnum] - horizstart[horizcellnum] + 1
            cell_nrows <- vertend[vertcellnum] - vertstart[vertcellnum] + 1
            indices <- matrix(1:(cell_nrows*cell_ncols), nrow=cell_nrows, 
                              ncol=cell_ncols)
            samp_indices <- sample(indices, nsamp)
            if (class(x) == 'RasterLayer') {
                xcell <- getValuesBlock(x, row=vertstart[vertcellnum],
                                        nrows=cell_nrows,
                                        col=horizstart[horizcellnum],
                                        ncols=cell_ncols)
            } else {
                xcell <- x[vertstart[vertcellnum]:vertend[vertcellnum],
                           horizstart[horizcellnum]:horizend[horizcellnum]]
            }
            retvals[retval_index:(retval_index + length(samp_indices) - 1)] <- xcell[samp_indices]
            retval_index <- retval_index + length(samp_indices)
        }
    }
    return(retvals)
}
