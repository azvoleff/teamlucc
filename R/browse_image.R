#' Simple function to output browse image as PNG
#'
#' @export
#' @param x input \code{RasterBrick} or \code{RasterStack} with at least three 
#' bands
#' @param filename a filename ending in .png
#' @param an optional mask \code{RasterLayer} to output below the browse image
#' @param maxpixels maximum number of pixels to use in plotting
#' @param DN_min minimum DN value
#' @param DN_max maximum DN value
#' @param r index in \code{x} of the band to use as the red band
#' @param g index in \code{x} of the band to use as the green band
#' @param b index in \code{x} of the band to use as the blue band
#' @param with width in pixels of output PNG file
#' @param height height in pixels of output PNG file
#' @param x_fun an optional function to apply to \code{x} after x is resampled 
#' according to \code{maxpixels}
#' @param m_fun an optional function to apply to \code{m} after m is resampled 
#' according to \code{maxpixels}
#' @return nothing - used for side-effect of saving browse image
browse_image <- function(x, filename, m=NULL, maxpixels=500000, DN_min=0, 
                         DN_max=255, r=3, g=2, b=1, width=1200, height=1200, 
                         x_fun=NULL, m_fun=NULL) {
    if (!tolower(extension(filename)) == '.png') {
        stop('output filename does not end in ".png"')
    }
    if (!file_test('-d', dirname(filename))) {
        stop('output path does not exist')
    }
    stopifnot(nlayers(m) == 1)

    x <- sampleRegular(x, size=maxpixels, asRaster=TRUE, useGDAL=TRUE)

    x <- calc(x, fun=function(vals) {
        vals[vals < DN_min] <- DN_min
        vals[vals > DN_max] <- DN_max
        vals <- ((vals - DN_min) / (DN_max - DN_min)) * 255
        if (!is.null(x_fun)) {
            vals <- x_fun(vals)
        }
        return(vals)
    }, datatype='INT1U')

    if (!is.null(m_fun)) {
        m <- calc(m, fun=m_fun, datatype=dataType(m))
    }

    png(filename=browse_file, width=width, height=height)
    if (!is.null(m)) {
        m <- sampleRegular(m, size=maxpixels, asRaster=TRUE, useGDAL=TRUE)
        if (nrow(x) > ncol(x)) par(mfrow=c(1, 2))
        else par(mfrow=c(2, 1))
        plotRGB(x, r=r, g=g, b=b, maxpixels=maxpixels)
        plot(m, maxpixels=maxpixels, axes=FALSE, legend=FALSE, box=FALSE)
    } else {
        plotRGB(x, r=r, g=g, b=b, maxpixels=maxpixels)
    }
    dev.off()
}
