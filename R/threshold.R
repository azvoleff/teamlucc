#' Automatically determine value for image thresholding
#'
#' The only method currently implemented is Huang's fuzzy thresholding method.  
#' The code for Huang's method was ported to C++ for \code{teamlucc} from the 
#' code in the Auto_threshold imageJ plugin by Gabriel Landini. See original 
#' code at:
#' http://www.mecourse.com/landinig/software/autothreshold/autothreshold.html
#'
#' This function will run in parallel if a parallel backend is registered with 
#' \code{\link{foreach}}.
#'
#' @export
#' @import foreach
#' @importFrom iterators iter
#' @param x the input image, as a matrix or raster
#' @param method the thresholding method. Currently only "huang" is 
#' implemented.
#' @param n_bin number of bins to use when calculating histogram
#' @param maxpixels maximum number of pixels size to use when calculating 
#' histogram
#' @return integer threshold value
#' @references Huang, L.-K., and M.-J. J. Wang. 1995. Image thresholding by 
#' minimizing the measures of fuzziness. Pattern recognition 28 (1):41--51.
threshold <- function(x, method="huang", n_bin=1000, maxpixels=5e5) {
    stopifnot(method %in% c("huang"))

    if (ncell(x) > maxpixels) {
        x <- sampleRegular(x, maxpixels, useGDAL=TRUE, asRaster=TRUE)
    }

    x <- stack(x)
    x <- setMinMax(x)
    mins <- minValue(x)
    maxs <- maxValue(x)

    bys <- (maxs - mins) / (n_bin)

    bandnum=minval=maxval=NULL
    thresholds <- foreach(bandnum=iter(1:nlayers(x)), minval=iter(mins), 
                          maxval=iter(maxs), by=iter(bys),
                          .packages=c('teamlucc'),
                          .combine=c) %dopar% {
        image_hist <- hist(x[[bandnum]], breaks=seq(minval, maxval+by, by=by), 
                           plot=FALSE, maxpixels=maxpixels)
        threshold_index <- threshold_Huang(image_hist$counts)
        image_hist$breaks[threshold_index]
    }

    return(thresholds)
}
