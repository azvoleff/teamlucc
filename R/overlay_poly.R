#' Make plot of image with overlaid polygon
#'
#' Useful for quick plots showing overlap between an area of interest (AOI) and 
#' a satellite image.
#'
#' @export
#' @import ggplot2
#' @importFrom rgeos gBuffer
#' @importFrom grid unit
#' @importFrom plyr join
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
overlay_poly <- function(x, y, out, title='', width=6.5, height=6.5, dpi=300, ...) {
    if (!(nlayers(x) %in% c(1, 3))) {
        stop('x must be a one or three band image')
    }

    y <- spTransform(y, CRS(proj4string(x)))
    x <- crop(x, y)
    x <- mask(x, y)

    # Resample the image based on the browse output DPI and image size. No need 
    # to stretch the entire Landsat image just to output a small browse image.
    agg_fact <- floor(max(nrow(x), ncol(x)) / (dpi * max(width, height)))
    if (agg_fact > 1) {
        new_params <- raster(x)
        extent(new_params) <- extent(x)
        res(new_params) <- res(x) * agg_fact
        x <- resample(x, new_params, method='ngb')
    }

    # Setup the zoi dataframe for ggplot2
    y@data$id <- rownames(y@data)
    y.points <- fortify(y, region="id")
    y.df <- join(y.points, y@data, by="id")

    # Setup color matrix for ggplot2, applying a linear 2 percent stretch
    x_vals <- linear_stretch(getValues(x))
    x_vals[is.na(x_vals)] <- 1
    if (nlayers(x) == 3) {
        fc_df <- data.frame(x=coordinates(x)[, 1],
                            y=coordinates(x)[, 2],
                            color_val=rgb(x_vals))
    } else if (nlayers(x) == 1) {
        fc_df <- data.frame(x=coordinates(x)[, 1],
                            y=coordinates(x)[, 2],
                            color_val=x_vals)
    } else {
        stop('overlay_poly only supports one or three band imagery')
    }

    # Fix for R CMD CHECK:
    color_val=long=lat=NULL

    theme_set(theme_bw(base_size=8))
    ggplot(fc_df) +
        geom_raster(aes(x, y, fill=color_val)) + coord_fixed() + 
        scale_fill_identity() + 
        labs(title=title) +
        theme(axis.text.x=element_blank(), axis.text.y=element_blank(),
              axis.title.x=element_blank(), axis.title.y=element_blank(),
              panel.background=element_blank(), panel.border=element_blank(),
              panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
              plot.background=element_blank(), axis.ticks=element_blank(),
              plot.margin=unit(c(.1, .1, .1, .1), 'cm')) +
        geom_path(data=y.df, aes(long, lat), color='blue', size=1)
    ggsave(out, width=width, height=height, dpi=dpi)

}
