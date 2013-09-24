#' Extract training data for use in a classification
#'
#' @export
#' @param x A \code{Raster*} object from which training data will be 
#' extracted. Training data will be extracted from each layer in a 
#' \code{RasterBrick} or \code{RasterStack}.
#' @param y A \code{SpatialPolygonsDataFrame} with training data polygons.
#' @return A data.frame with the training data. Each row will contain the 
#' response ('y') as the first column, with the remaining columns containing 
#' the values at that location of each band in the raster stack. 
#' @examples
#' x <- stack(system.file('extdata/L5TSR_1986.dat', package='teamr'))
#' data(L5TSR_1986_training)
#' train_data <- extract_training_data(x, L5TSR_1986_training)
extract_training_data <- function(x, y) {
    if (projection(x) != projection(y)) {
        stop('Coordinate systems do not match')
    }
    pixels <- extract(x, y, small=TRUE, df=TRUE)
    # Add column with class labels
    y$ID <- seq(1, nrow(y))
    pixels <- cbind(y=y$Poly_Type[match(pixels$ID, y$ID)], 
                    pixels)
    pixels <- pixels[!(names(pixels) == 'ID')]
    return(pixels)
}
