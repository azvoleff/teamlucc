#' Apply a raster function with edge effects over a series of blocks
#'
#' This function can be useful when applying windowed functions over a raster, 
#' as with \code{glcm}. This function allows windows functions that have edge 
#' effects to be applied over a raster in block-by-block fashion.  
#' \code{applyWindowed} avoids the striping that would result if the edge 
#' effects were ignored.
#'
#' @export
#' @param x a \code{Raster*}
#' @param fun the function to apply
#' @param edge length 2 numberic with number of rows on top and bottom with 
#' edge effects, defined as c(top, bottom)
#' @param chunksize the number of rows to read per block (passed to 
#' \code{raster} \code{blockSize} function.
#' @param ... additional arguments to pass to \code{fun}
#' @examples
#' \dontrun{
#' L5TSR_1986_b1 <- raster(L5TSR_1986, layer=1)
#' min_x <- cellStats(L5TSR_1986_b1, 'min')
#' max_x <- cellStats(L5TSR_1986_b1, 'max')
#' applyWindowed(L5TSR_1986_b1, glcm, edge=c(1, 3), min_x=min_x, max_x=max_x)
#' }
applyWindowed <- function(x, fun, edge=c(0, 0), chunksize=NULL, ...) {
    if ((length(edge) != 2) || (class(edge) != 'numeric') || any(edge < 0)) {
        stop('edge must be a length 2 positive numeric')
    }

    if (is.null(chunksize)) {
        bs <- blockSize(x)
    } else {
        bs <- blockSize(x, chunksize)
    }
    n_blocks <- bs$n

    # bs_mod is the blocksize that will contain blocks that have been expanded 
    # to avoid edge effects
    bs_mod <- bs
    # Expand blocks to account for edge effects on the top:
    bs_mod$row[2:n_blocks] <- bs_mod$row[2:n_blocks] - edge[1]
    # Need to read additional rows from these blocks to avoid an offset
    bs_mod$nrows[2:n_blocks] <- bs_mod$nrows[2:n_blocks] + edge[1]

    # Read additional bottom rows to account for edge effects on the bottom:
    bs_mod$nrows[1:(n_blocks - 1)] <- bs_mod$nrows[1:(n_blocks - 1)] + edge[2]

    if (any(bs_mod$row < 1)) {
        stop('too many blocks to read without edge effects - try increasing chunksize')
    } else if (any((bs_mod$nrows + bs_mod$row - 1) > nrow(x))) {
        stop('too many blocks to read without edge effects - try increasing chunksize')
    }
    
    started_writes <- FALSE
    for (block_num in 1:bs$n) {
        this_block <- getValues(x, row=bs_mod$row[block_num], nrows=bs_mod$nrows[block_num],
                                format='matrix')
        out_block <- fun(this_block, ...)
        layer_names <- dimnames(out_block)[[3]]
        # Drop the padding added to top to avoid edge effects, unless we are 
        # really on the top of the image, where top edge effects cannot be 
        # avoided
        if ((block_num != 1) && (edge[1] > 0)) {
            out_block <- out_block[-(1:edge[1]), , ]
            # The below line is needed to maintain a 3 dimensional array, 
            # even when an n x m x 1 array is returned from 
            # calc_texture_full_image because a single statistic was chosen. 
            # Without the below line, removing a row will coerce the 3d array 
            # to a 2d matrix, and the bottom padding removal will fail as it 
            # references a 3d matrix).
            if (length(dim(out_block)) < 3) dim(out_block) <- c(dim(out_block), 1)
        }
        # Drop the padding added to bottom to avoid edge effects, unless we are 
        # really on the bottom of the image, where bottom edge effects cannot 
        # be avoided
        if ((block_num != n_blocks) && (edge[2] > 0)) {
            out_block <- out_block[-((nrow(out_block)-edge[2]+1):nrow(out_block)), , ]
            if (length(dim(out_block)) < 3) dim(out_block) <- c(dim(out_block), 1)
        }
        if (!started_writes) {
            # Setup an output raster with number of layers equal to the number 
            # of layers in out_block, and extent/resolution equal to extent and 
            # resolution of x
            if (dim(out_block)[3] == 1) {
                out <- raster(x)
            } else {
                out <- brick(stack(rep(c(x), dim(out_block)[3])), values=FALSE)
            }
            out <- writeStart(out, filename=rasterTmpFile())
            names(out) <- layer_names
            started_writes <- TRUE
        }
        # To write to a RasterBrick the out_block needs to be structured as 
        # a 2-d matrix with bands in columns and columns as row-major vectors
        if (dim(out_block)[3] == 1) {
            out_block <- aperm(out_block, c(3, 2, 1))
            out_block <- matrix(out_block, ncol=nrow(out_block))
        } else {
            out_block <- aperm(out_block, c(3, 2, 1))
            out_block <- matrix(out_block, ncol=nrow(out_block), byrow=TRUE)
        }
        out <- writeValues(out, out_block, bs$row[block_num])
    }
    out <- writeStop(out)
    return(out)
}
