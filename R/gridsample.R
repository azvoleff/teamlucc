#' Draw a random sample from a grid laid out on a RasterLayer or matrix
#'
#' This function is used to subsample a \code{RasterLayer} or \code{matrix} by 
#' dividing the dataset into a grid of \code{horizcells} x \code{vertcells}, 
#' and by then drawing a sample of size \code{nsamp} from within each grid 
#' cell.
#'
#' @export
#' @param x a matrix or RasterLayer to draw sample from
#' @param horizcells how many cells to break the raster in horizontally (over 
#' the columns)
#' @param vertcells how many cells to break the raster in vertically (over 
#' the rows)
#' @param nsamp how many samples to draw from each grid cell
#' @param rowmajor whether to return indices in row-major format (default is 
#' column-major). Row-major format is useful in conjunction with \code{Raster*} 
#' objects.
#' @param replace whether to sample with replacement (within each grid cell)
#' @return vector of sample indices
#' @note TODO: Recode in C++ for speed.
#' @examples
#' # Make a 100x100 matrix
#' x <- matrix(1:10000, nrow=100)
#' # Subsample the matrix by breaking it up into a 10 x 10 grid, and take 10
#' # random samples from each grid cell without replacement (these are the
#' # default parameters).
#' y <- gridsample(x)
gridsample <- function(x, horizcells=10, vertcells=10, nsamp=10, 
                       rowmajor=FALSE, replace=FALSE) {
    # horizstart is a vector of column numbers of the first column in each cell 
    # in the grid
    horizstart <- round(seq(1, ncol(x), ncol(x) / horizcells))
    # horizend is a vector of column numbers of the last column in each cell in 
    # the grid
    if (length(horizstart) > 1) {
        horizend <- c((horizstart - 1)[2:length(horizstart)], ncol(x))
    } else {
        horizend <- c(ncol(x))
    }
    # vertstart is a vector of row numbers of the first row in each cell in the 
    # grid
    vertstart <- round(seq(1, nrow(x), nrow(x) / vertcells))
    # vertend is a vector of row numbers of the last row in each cell in the
    # grid
    if (length(vertstart) > 1) {
        vertend <- c((vertstart - 1)[2:length(vertstart)], nrow(x))
    } else {
        vertend <- c(nrow(x))
    }
    # retvalues is a vector to store the sample values chosen from x
    # sampindices is a vector to store the column major indices of the locations 
    # of each sampled value in x
    sampindices <- vector('numeric', horizcells * vertcells * nsamp)
    # retval_index tracks position in the vector storing the sampled values 
    # returned from this function
    retval_index <- 1
    # cell1row is used in calculating the column-major index of the first entry 
    # (top left) of each grid cell
    cell1row <- 1
    for (vertcellnum in 1:length(vertstart)) {
        # cell1col is used in calculating the column-major index of the first 
        # entry (top left) of each grid cell
        cell1col <- 1
        # cell_nrows is the number of rows in this particular cell (cells may 
        # have varying numbers of rows due to rounding)
        cell_nrows <- vertend[vertcellnum] - vertstart[vertcellnum] + 1
        for (horizcellnum in 1:length(horizstart)) {
            # cell1colmajindex is the column-major index of the first (top 
            # left) value in this particular grid cell
            cell1colmajindex <- cell1row + nrow(x) * (cell1col - 1)
            # cell_ncols is the number of rows in this particular cell (cells 
            # may have varying numbers of columns due to rounding)
            cell_ncols <- horizend[horizcellnum] - horizstart[horizcellnum] + 1
            # cell_colmaj_indices is a matrix of column-major indices of the 
            # position of each value in this grid cell, within the larger 
            # matrix x
            cell_colmaj_indices <- matrix(rep(1:cell_nrows, cell_ncols),
                                          nrow=cell_nrows) +
                                   matrix(rep(seq(cell1colmajindex, by=nrow(x), 
                                                  length.out=cell_ncols), 
                                              cell_nrows), nrow=cell_nrows, 
                                          byrow=TRUE) - 1
            samp_indices <- sample(cell_colmaj_indices, nsamp, replace=replace)
            sampindices[retval_index:(retval_index + length(samp_indices) - 1)] <- samp_indices
            retval_index <- retval_index + length(samp_indices)
            cell1col  <- cell1col + cell_ncols
        }
        cell1row <- cell1row + cell_nrows
    }
    if (rowmajor) {
        sampindices <- ((sampindices - 1) %% nrow(x)) * ncol(x) + ((sampindices - 1) %/% nrow(x) + 1)
    }
    return(sampindices)
}
