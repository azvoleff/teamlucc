#' Exports statistics on pixels within each of a set of land cover classes
#'
#' @export
#' @import plyr
#' @importFrom reshape2 melt
#' @param x A \code{RasterLayer} from which class statistics will be 
#' calculated.
#' @param y A \code{SpatialPolygonsDataFrame} with cover class 
#' polygons
#' @param classcol the name of the column containing the response variable (for 
#' example the land cover type of each pixel)
#' @return A data.frame of class statistics.
#' @examples
#' L5TSR_1986 <- stack(system.file('extdata/L5TSR_1986.dat', package='teamr'))
#' class_statistics(L5TSR_1986, L5TSR_1986_2001_training, "t1_class")
class_statistics <- function(x, y, classcol) {
    if (projection(x) != projection(y)) {
        stop('Coordinate systems do not match')
    }
    if (class(y) == "SpatialPolygonsDataFrame") {
        pixels <- extract_training_data(x, y, classcol)
    } else if (class(y) %in% c("RasterLayer", "RasterBrick", 
                                         "RasterStack")) {
        stop('Error: class_statistics cannot yet handle Raster* objects')
    }
    pixels <- melt(pixels, idvar='y')
    # Set y and variable to NULL to pass R CMD CHECK without notes
    value=variable=NULL
    class_stats <- ddply(pixels, .(y, variable), summarize,
                         r_mean=mean(value),
                         r_sd=sd(value),
                         r_min=min(value),
                         r_max=max(value),
                         n_pixels=length(value))
    return(class_stats)
}
