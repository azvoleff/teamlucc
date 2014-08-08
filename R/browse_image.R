plotprep <- function(x, maxpixels=500000, DN_min=0, DN_max=255, x_fun=NULL) {
    if (ncell(x) > maxpixels) {
        x <- sampleRegular(x, size=maxpixels, asRaster=TRUE, useGDAL=TRUE)
    }
    x <- calc(x, fun=function(vals) {
        vals[vals < DN_min] <- DN_min
        vals[vals > DN_max] <- DN_max
        vals <- ((vals - DN_min) / (DN_max - DN_min)) * 255
        if (!is.null(x_fun)) {
            vals <- x_fun(vals)
        }
        return(vals)
    }, datatype='INT1U')
    return(x)
}

#' Simple function to make small preview plot from large raster image
#'
#' @export
#' @param x input \code{RasterBrick} or \code{RasterStack} with at least three 
#' bands
#' @param m an optional mask \code{RasterLayer} to output below the browse 
#' image
#' @param maxpixels maximum number of pixels to use in plotting
#' @param DN_min minimum DN value
#' @param DN_max maximum DN value
#' @param r index in \code{x} of the band to use as the red band
#' @param g index in \code{x} of the band to use as the green band
#' @param b index in \code{x} of the band to use as the blue band
#' @param x_fun an optional function to apply to \code{x} after x is resampled 
#' according to \code{maxpixels}
#' @param m_fun an optional function to apply to \code{m} after m is resampled 
#' according to \code{maxpixels}
#' @return nothing - used for side-effect of saving browse image
browse_image <- function(x, m=NULL, maxpixels=500000, DN_min=0, DN_max=255, 
                         r=3, g=2, b=1, x_fun=NULL, m_fun=NULL) {
    if (!is.null(m)) stopifnot(nlayers(m) == 1)

    x <- plotprep(x, maxpixels=500000, DN_min=DN_min, DN_max=DN_max, 
                  x_fun=x_fun)

    if (!is.null(m) && !is.null(m_fun)) {
        m <- calc(m, fun=m_fun, datatype=dataType(m))
    }

    if (!is.null(m)) {
        m <- sampleRegular(m, size=maxpixels, asRaster=TRUE, useGDAL=TRUE)
        if (nrow(x) > ncol(x)) par(mfrow=c(1, 2))
        else par(mfrow=c(2, 1))
        plotRGB(x, r=r, g=g, b=b, maxpixels=maxpixels)
        plot(m, maxpixels=maxpixels, axes=FALSE, legend=FALSE, box=FALSE)
    } else {
        plotRGB(x, r=r, g=g, b=b, maxpixels=maxpixels)
    }
}
