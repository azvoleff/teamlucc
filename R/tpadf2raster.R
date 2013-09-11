#' Function to convert TIMESAT .tpa data.frame to an R raster. 
#'
#' @export
#' @param tpadf A TPA data.frame as output by tpa2df
#' @param base_image_file A string giving the location of a raster file to use 
#' for georeferencing the output raster. Use one of the original raster files 
#' that was input to TIMESAT.
#' @param variable A string giving the variable name to write to a raster.  Can 
#' be one of: start, end, length, base_value, peak_time, peak_value, amp, 
#' left_deriv, right_deriv, large_integ, and small_integ.
#' @return A raster object
#' @examples
#' # TODO: Need to add examples here, and need to include a sample TIMESAT tpa 
#' # file in the package data.
tpadf2raster <- function(tpadf, base_image_file, variable) {
    # Parameters:
    # base_image_file should be one of the original MODIS files that was fed 
    # into TIMESAT.
    #
    # variable can be any of the variables in the tpa dataframe.
    require(raster)

    if (missing(tpadf) || !is.data.frame(tpadf)) {
        stop('must specify a tpa data.frame')
    } else if (missing(base_image_file) || !file.exists(base_image_file)) {
        stop('must specify a valid base image raster')
    }

    var_col <- grep(paste('^', variable, '$', sep=''), names(tpadf))
    if (length(var_col) == 0) {
        stop(paste(variable, 'not found in tpa dataframe'))
    }
    base_image <- raster(base_image_file)
    ncol(base_image) * nrow(base_image) * 2

    seasons <- sort(unique(tpadf$season))
    out_rasters <- c()
    for (season in sort(unique(tpadf$season))) {
        season_data <- tpadf[tpadf$season == season, ]
        data_matrix <- matrix(NA, nrow(base_image), ncol(base_image))
        vector_indices <- (nrow(data_matrix) * season_data$col) - 
            (nrow(data_matrix) - season_data$row)
        data_matrix[vector_indices] <- season_data[, var_col]
        out_raster <- raster(data_matrix, template=base_image)
        out_rasters <- c(out_rasters, out_raster)
    }
    out_rasters <- stack(out_rasters)
    names(out_rasters) <- seasons

    return(out_rasters)
}
