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
#' @param out_file path to save output raster to
#' @param edge length 2 numberic with number of rows on top and bottom with 
#' edge effects, defined as c(top, bottom)
#' @param ... additional arguments to pass to \code{fun}
#' @examples
#' \dontrun{
#' applyWindowed(L5TSR_1986, glcm, 'out.envi', edge=c(1, 3))
#' }
applyWindowed <- function(x, fun, out_file, edge=c(0, 0), ...) {
    if ((length(edge) != 2) || (class(edge) != 'numeric') || any(edge < 0)) {
        stop('edge must be a length 2 positive numeric')
    }

    # bs is the blocksize that will track writes
    bs <- blockSize(x)
    n_blocks <- bs$n

    # bs_mod is the blocksize that will contain blocks that have been expanded 
    # to avoid edge effects
    bs_mod <- bs
    # Expand blocks to account for edge effects on the top:
    bs_mod$row[2:n_blocks] <- bs_mod$row[2:n_blocks] - edge[1]
    # Need to read an additional row from these blocks to avoid an offset
    bs_mod$nrows[2:n_blocks] <- bs_mod$nrows[2:n_blocks] + edge[1]

    # Read additional bottom rows to account for edge effects on the bottom:
    bs_mod$nrows[1:(n_blocks - 1)] <- bs_mod$nrows[1:(n_blocks - 1)] + edge[2]
    
    started_writes <- FALSE
    for (block_num in 1:bs$n) {
        this_block <- getValues(x, row=bs_mod$row[block_num], nrows=bs_mod$nrows[block_num])
        out_block <- fun(this_block, ...)
        # Drop the extra edge rows added to avoid edge effects:
        if (edge[1] > 0) out_block <- out_block[-edge[1]]
        if (edge[2] > 0) out_block <- out_block[-((nrow(out_block)-edge[1]):nrow(out_block))]
        if (!started_writes) {
            # Setup an output raster with number of layers equal to the number 
            # of layers in out_block, and extent/resolution equal to extent and 
            # resolution of x
            if (length(dim(out_block)) <= 2) {
                out <- raster(x)
            } else {
                out <- stack(rep(c(x), dim(out_block)[3]))
            }
            out <- writeStart(out, out_file)
            started_writes <- TRUE
        }
        writeValues(out, out_block, bs$row[block_num])
    }
    out <- writeStop(out)

    return(out)
}
