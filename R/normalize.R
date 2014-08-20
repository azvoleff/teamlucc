#' Normalizes two rasters
#'
#' Performs relative normalization on two rasters using model II regression. 
#' Based on the approach in the \code{relnorm} function in the \code{landsat} 
#' package.
#'
#' This function will run in parallel if a parallel backend is registered with 
#' \code{\link{foreach}}.
#'
#' @export
#' @importFrom iterators iter
#' @import foreach
#' @importFrom lmodel2 lmodel2
#' @param x a \code{Raster*} to use as the base image
#' @param y a \code{Raster*} to normalize to the base image
#' @param msk a \code{RasterLayer} with missing values in \code{x} or in {y} 
#' coded as 1, and all other values coded as 0 (optional)
#' @param method the regression method to use (must be a method recognized by 
#' \code{lmodel2}
#' @param size the number of pixels to use in developing the model
#' @return a \code{Raster*} of \code{y} normalized to \code{x}
#' @examples
#' L5TSR_2001_normed_1 <- normalize(L5TSR_1986, L5TSR_2001)
#' plotRGB(L5TSR_2001_normed_1, stretch='lin')
#'
#' # Use only half as many pixels to calculate the models
#' L5TSR_2001_normed_2 <- normalize(L5TSR_1986, L5TSR_2001, 
#'                                  size=ncell(L5TSR_1986)/2)
#' plotRGB(L5TSR_2001_normed_2, stretch='lin')
#' @references
#' Sarah Goslee. Analyzing Remote Sensing Data in {R}: The {landsat} Package.  
#' Journal of Statistical Software, 2011, 43:4, pg 1--25.  
#' http://www.jstatsoft.org/v43/i04/
normalize <- function(x, y, msk, method="MA", size=ncell(x)) {
    orig_datatype <- dataType(y)[1]
    compareRaster(x, y)
    stopifnot(nlayers(x) == nlayers(y))
    stopifnot(size <= ncell(x))
    if (!missing(msk)) {
        compareRaster(x, msk)
        stopifnot(nlayers(msk) == 1)
    }

    if (size < ncell(x)) {
        # Note that sampleRegular with cells=TRUE returns cell numbers in the 
        # first column
        x_vals <- sampleRegular(x, size=size, cells=TRUE)
        if (!missing(msk)) {
            x_vals <- x_vals[!(msk[x_vals[, 1]]), ]
        }
        y_vals <- y[x_vals[, 1]]
        x_vals <- x_vals[, -1]
    } else {
        x_vals <- getValues(x)
        y_vals <- getValues(y)
        if (!missing(msk)) {
            x_vals <- x_vals[!getValues(msk), ]
            y_vals <- y_vals[!getValues(msk), ]
        }
    }

    names(y_vals) <- names(x_vals)

    if (nlayers(y) > 1) {
        unnormed_layer <- x_sample <- y_sample <- NULL
        normed_y <- foreach(unnormed_layer=unstack(y),
                            x_sample=iter(x_vals, by='column'),
                            y_sample=iter(y_vals, by='column'),
                            .combine='addLayer', .multicombine=TRUE, 
                            .init=raster(),
                            .packages=c('raster', 'lmodel2', 'rgdal')) %dopar% {
            model <- suppressMessages(lmodel2(x_sample ~ y_sample, nperm=0))
            model <- model$regression.results[model$regression.results[, "Method"] == method, ]
            names(model) <- gsub("^ *", "", names(model))
            normed_layer <- model$Slope * unnormed_layer + model$Intercept
        }
    } else {
        model <- suppressMessages(lmodel2(x_vals ~ y_vals, nperm=0))
        model <- model$regression.results[model$regression.results[, "Method"] == method, ]
        names(model) <- gsub("^ *", "", names(model))
        normed_y <- model$Slope * y + model$Intercept
    }

    if (!missing(msk)) {
        # Copy masked values back into the output raster
        normed_y[msk] <- y[msk]
    }

    #dataType(normed_y) <- orig_datatype

    return(normed_y)
}
