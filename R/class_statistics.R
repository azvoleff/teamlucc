#' Exports statistics on pixels within each of a set of land cover classes
#'
#' @export
#' @param rast A \code{RasterLayer} from which class statistics will be 
#' calculated.
#' @param shp A shapefile with class polygons
#' @return A data.frame of class statistics.
#' @examples
#' # TODO: Add example
#' # library(teamr)
#' # shp <- 'H:/Data/TEAM/VB/Vectors/VB_training_1986_037_LT5.shp'
#' # img_file <- 'H:/Data/TEAM/VB/Rasters/Landsat/1986_037_LT5/proc/lndsr.LT50150531986037XXX18.bsq'
#' # rast <- raster(img_file, band=4)
#' # class_statistics(rast, shp)
class_statistics <- function(rast, shp) {
    layer_name <- gsub('.shp$', '', basename(shp), ignore.case=TRUE)
    class_polys <- readOGR(dirname(shp), layer_name)
    if (projection(rast) != projection(class_polys)) {
        stop('Coordinate systems do not match')
    }
    class_stats <- data.frame(Poly_Type=unique(class_polys$Poly_Type))
    class_stats$r_mean <- NA
    class_stats$r_sd <- NA
    class_stats$r_min <- NA
    class_stats$r_max <- NA
    class_stats$n_polys <- NA
    class_stats$mean_n_per_poly <- NA
    class_stats$sd_n_per_poly <- NA
    class_stats$min_n_per_poly <- NA
    class_stats$max_n_per_poly <- NA
    for (n in 1:nrow(class_stats)) {
        these_polys <- class_polys[class_polys$Poly_Type == class_stats$Poly_Type[n], ]
        pixels <- extract(rast, these_polys, small=TRUE)
        class_stats$r_mean[n] <- mean(unlist(pixels))
        class_stats$r_sd[n] <- sd(unlist(pixels))
        class_stats$r_min[n] <- min(unlist(pixels))
        class_stats$r_max[n] <- max(unlist(pixels))
        class_stats$n_polys[n] <- length(these_polys)
        class_stats$mean_n_per_poly[n] <- mean(sapply(pixels, function(x) length(!is.na(x))))
        class_stats$sd_n_per_poly[n] <- sd(sapply(pixels, function(x) length(!is.na(x))))
        class_stats$min_n_per_poly[n] <- min(sapply(pixels, function(x) length(!is.na(x))))
        class_stats$max_n_per_poly[n] <- max(sapply(pixels, function(x) length(!is.na(x))))
    }
    return(class_stats)
}
