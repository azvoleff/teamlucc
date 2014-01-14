.do_stretch <- function(x, pct, max_val) {
    lower <- quantile(x, prob=0 + pct/100, na.rm=TRUE)
    upper <- quantile(x, prob=1 - pct/100, na.rm=TRUE)
    x[x > upper] <- upper
    x[x < lower] <- lower
    x <- ((x - lower) / (upper-lower)) * max_val
    return(x)
}

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
linear_stretch <- function(x, pct=2, max_val=1) {
    # Applies linear stretch (2 percent by default). Assumes image is arranged 
    # with bands in columns. Returns the image with stretch applied and bands 
    # rescaled to range from 0 - max_val.
    if ((pct < 0) | pct >= 50) {
        stop('pct must be > 0 and < 50')
    }
    if (class(x) %in% c('RasterLayer')) {
        x <- setValues(x, .do_stretch(getValues(x), pct, max_val))
        return(x)
    } else if (class(x) %in% c('RasterStack', 'RasterBrick')) {
        for (n in 1:nlayers(x)) {
            x <- setValues(x, .do_stretch(getValues(raster(x, layer=n)),
                                          pct, max_val), layer=n)
        }
        return(x)
    } else if (is.null(dim(x)) || (length(dim(x)) == 2) || (dim(x)[3] == 1)) {
        # Handle a vector, matrix, or n x m x 1 array
        return(.do_stretch(x, pct, max_val))
    }
}
