#' Function to convert TIMESAT tts data.frame an R raster.
#'
#' @export
#' @param x A TTS data.frame as output by tts2df
#' @param base_image A string giving the location of a raster file to use 
#' for georeferencing the output raster. Use one of the original raster files 
#' that was input to TIMESAT.
#' @return A raster object
#' @examples
#' # TODO: Need to add examples here, and need to include a sample TIMESAT tpa 
#' # file in the package data.
ttsdf2raster <- function(x, base_image) {
    if (missing(x) || !is.data.frame(x)) {
        stop('must specify a tts data.frame')
    } else if (missing(base_image) || !file.exists(base_image)) {
        stop('must specify a valid base image raster')
    }

    t_cols <- grep('^t[0-9]{1,4}$', names(x))

    base_image <- raster(base_image)

    out_rasters <- c()
    for (t_col in t_cols) {
        this_time_data <- x[, t_col]
        data_matrix <- matrix(NA, nrow(base_image), ncol(base_image))
        vector_indices <- (nrow(data_matrix) * x$col) - 
            (nrow(data_matrix) - x$row)
        data_matrix[vector_indices] <- this_time_data
        out_raster <- raster(data_matrix, template=base_image)
        out_rasters <- c(out_rasters, out_raster)
    }
    out_rasters <- stack(out_rasters)
    names(out_rasters) <- names(x[t_cols])

    return(out_rasters)
}
