#' Exports statistics on pixels within each of a set of land cover classes
#'
#' @export
#' @importFrom dplyr group_by summarize
#' @importFrom reshape2 melt
#' @param x A \code{RasterLayer} from which class statistics will be 
#' calculated.
#' @param y A \code{SpatialPolygonsDataFrame} with cover class 
#' polygons
#' @param class_col the name of the column containing the response variable 
#' (for example the land cover type of each pixel)
#' @return A data.frame of class statistics.
#' @examples
#' class_statistics(L5TSR_1986, L5TSR_1986_2001_training, "class_1986")
class_statistics <- function(x, y, class_col) {
    if (projection(x) != projection(y)) {
        stop('Coordinate systems do not match')
    }
    if (class(y) == "SpatialPolygonsDataFrame") {
        pixels <- get_pixels(x, y, class_col)
    } else if (class(y) %in% c("RasterLayer", "RasterBrick", 
                                         "RasterStack")) {
        stop('class_statistics cannot yet handle Raster* objects')
    }
    pixels <- melt(data.frame(pixels@x, y=pixels@y), idvar='y')
    # Set y and variable to NULL to pass R CMD CHECK without notes
    value=variable=NULL
    class_stats <- summarize(group_by(pixels, y, variable), mean=mean(value), 
                             sd=sd(value), min=min(value), max=max(value), 
                             n_pixels=length(value))
    class_stats <- class_stats[order(class_stats$variable, class_stats$y), ]
    return(class_stats)
}
