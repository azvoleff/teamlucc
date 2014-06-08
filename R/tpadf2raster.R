#' Function to convert TIMESAT .tpa data.frame to an R raster. 
#'
#' @export
#' @param x A TPA data.frame as output by tpa2df
#' @param base_image A string giving the location of a raster file to use 
#' for georeferencing the output raster. Use one of the original raster files 
#' that was input to TIMESAT.
#' @param variable A string giving the variable name to write to a raster.  Can 
#' be one of: start, end, length, base_value, peak_time, peak_value, amp, 
#' left_deriv, right_deriv, large_integ, and small_integ.
#' @return A raster object
#' @examples
#' # TODO: Need to add examples here, and need to include a sample TIMESAT tpa 
#' # file in the package data.
tpadf2raster <- function(x, base_image, variable) {
    if (missing(x) || !is.data.frame(x)) {
        stop('must specify a tpa data.frame')
    } else if (missing(base_image) || !file.exists(base_image)) {
        stop('must specify a valid base image raster')
    }

    var_col <- grep(paste('^', variable, '$', sep=''), names(x))
    if (length(var_col) == 0) {
        stop(paste(variable, 'not found in tpa dataframe'))
    }
    base_image <- raster(base_image)
    ncol(base_image) * nrow(base_image) * 2

    seasons <- sort(unique(x$season))
    out_rasters <- c()
    for (season in sort(unique(x$season))) {
        season_data <- x[x$season == season, ]
        data_matrix <- matrix(NA, nrow(base_image), ncol(base_image))
        vector_indices <- (nrow(data_matrix) * season_data$col) - 
            (nrow(data_matrix) - season_data$row)
        data_matrix[vector_indices] <- season_data[, var_col]
        out_raster <- raster(data_matrix, template=base_image)
        out_rasters <- c(out_rasters, out_raster)
    }
    out_rasters <- stack(out_rasters)
    names(out_rasters) <- paste0('season_', seasons)

    return(out_rasters)
}
