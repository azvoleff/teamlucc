#' Make plot of image with overlaid polygon
#'
#' Useful for quick plots showing overlap between an area of interest (AOI) and 
#' a satellite image.
#'
#' @export
#' @importFrom grid unit
#' @import ggplot2
#' @param x image as a \code{Raster*} object
#' @param y polygon to overlay, as \code{SpatialPolygonDataFrame}
#' @param out filename for output image. The extension of this file will 
#' determine the output file format (png, pdf, etc.).
#' @param title title of plot
#' @param width width (in inches) of output image
#' @param height height (in inches) of output image
#' @param dpi DPI for output image
#' @param ... additional arguments to \code{ggsave}
#' @return used for writing plot to a file
overlay_poly <- function(x, y, out, title='', width=4, height=4, dpi=150,
                         ...) {
    stop('overlay_poly is not yet supported')
    if (!(nlayers(x) %in% c(1, 3))) {
        stop('x must be a one or three band image')
    }

    y <- spTransform(y, CRS(proj4string(x)))

    maxpixels <- width*height*dpi^2
    # Below is based on gplot from rasterVis package
    nl <- nlayers(x)
    if (ncell(x)/nl > maxpixels) {
        x <- sampleRegular(x, maxpixels*nl, ext=y, asRaster=TRUE)
    }

    x <- mask(x, y)

    coords <- xyFromCell(x, seq_len(ncell(x)))
    dat <- data.frame(linear_stretch(getValues(x)))
    dat <- dat[complete.cases(dat), ]
    coords <- coords[complete.cases(dat), ]
    if (nlayers(x) == 3) {
        dat <- rgb(dat)
    }
    dat <- data.frame(coords, color=dat)

    color=long=lat=NULL  # fix for R CMD CHECK
    theme_set(theme_bw(base_size=8))
    p <- ggplot(dat) +
        geom_raster(aes(x, y, fill=color)) + coord_fixed() +
        scale_fill_identity() +
        labs(title=title) +
        theme(axis.text.x=element_blank(), axis.text.y=element_blank(),
              axis.title.x=element_blank(), axis.title.y=element_blank(),
              panel.background=element_blank(), panel.border=element_blank(),
              panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
              plot.background=element_blank(), axis.ticks=element_blank(),
              plot.margin=unit(c(.1, .1, .1, .1), 'cm'))
    p
    ggsave(out, width=width, height=height, dpi=dpi, ...)
}
