#' Apply a linear stretch to an image
#'
#' Applies a linear stretch to an image (default linear 2% stretch), and 
#' returns the image with each band individually stretched and rescaled to 
#' range between zero and \code{max_val} (default of \code{max_val} is 1).
#'
#' @export
#' @param x image to stretch
#' @param pct percent stretch
#' @param max_val maximum value of final output (image will be rescaled to 
#' range from 0 - \code{max_val})
#' @return image with stretch applied
linear_stretch <- function(x, pct=2, max_val=1, ...) {
    # Applies linear stretch (2 percent by default). Assumes image is arranged 
    # with bands in columns. Returns the image with stretch applied and bands 
    # rescaled to range from 0 - max_val.
    if ((pct < 0) | pct > 100) {
        stop('pct must be > 0 and < 100')
    }
    for (n in 1:ncol(x)) {
        pct2 <- quantile(x[, n], prob=0 + pct/100, na.rm=TRUE)
        pct98 <- quantile(x[, n], prob=1 - pct/100, na.rm=TRUE)
        x[, n] <- ((x[, n] - pct2) / pct98) * max_val
        x[, n][x[, n] < 0] <- 0
        x[, n][x[, n] > 1] <- 1
    }
    return(x)
}
