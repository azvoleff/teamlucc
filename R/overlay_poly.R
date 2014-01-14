#' Make plot of image with overlaid polygon
#'
#' Useful for quick plots showing overlap between an area of interest (AOI) and 
#' a satellite image.
#'
#' @export
#' @import ggplot2
#' @importFrom rgeos gBuffer
#' @importFrom plyr unit
#' @param x image as a \code{Raster*} object
#' @param out filename for output image. The extension of this file will 
#' determine the output file format (png, pdf, etc.).
#' @param width width (in inches) of output image
#' @param height height (in inches) of output image
#' @param dpi DPI for output image
#' @param ... additional arguments to \code{ggsave}
#' @return used for writing plot to a file
overlay_poly <- function(x, y, out, width=6.5, height=6.5, dpi=300, ...) {
    if (!(nlayers(x) %in% c(1, 3))) {
        stop('x must be a one or three band image')
    }

    # Resample the image based on the browse output DPI and image size. No need 
    # to stretch the entire Landsat image just to output a small browse image.
    agg_fact <- floor(max(nrow(x), ncol(x)) / (DPI * max(PLOT_WIDTH, 
                                                         PLOT_HEIGHT)))
    x <- aggregate(x, fact=agg_fact)
    fc_zoi <- spTransform(zoi, CRS(proj4string(x)))
    x <- crop(x, fc_zoi)
    x <- mask(x, fc_zoi)

    year <- image_list[n, ]$year
    julian_day <- image_list[n, ]$julian_day
    sensor <- image_list[n, ]$sensor
    main_title <- paste(year, '-', sprintf("%03i", julian_day),
                        ' (', sensor, ')', sep='')
    sub_title <- as.character(image_list[n, ]$file)
    full_title <- substitute(atop(main_title, atop(sub_title)), 
                             list(main_title=main_title, sub_title=sub_title))

    # Setup the zoi dataframe for ggplot2
    fc_zoi@data$id <- rownames(fc_zoi@data)
    fc_zoi.points <- fortify(fc_zoi, region="id")
    fc_zoi.df <- join(fc_zoi.points, fc_zoi@data, by="id")

    # Setup color matrix for ggplot2, applying a linear 2 percent stretch
    x_vals <- getValues(x)
    x_vals <- linear_stretch(x_vals)
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

    theme_set(theme_bw(base_size=8))
    ggplot(fc_df) +
        geom_raster(aes(x, y, fill=color_val)) + coord_fixed() + 
        scale_fill_identity() + 
        labs(title=full_title) +
        theme(axis.text.x=element_blank(), axis.text.y=element_blank(),
              axis.title.x=element_blank(), axis.title.y=element_blank(),
              panel.background=element_blank(), panel.border=element_blank(),
              panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
              plot.background=element_blank(), axis.ticks=element_blank(),
              plot.margin=unit(c(.1, .1, .1, .1), 'cm')) +
        geom_path(data=fc_zoi.df, aes(long, lat), color='blue', size=1)
    ggsave(out, width=width, height=height, dpi=dpi)
}
