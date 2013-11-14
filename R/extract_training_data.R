#' Extract training data for use in a classification
#'
#' @export
#' @param x a \code{Raster*} object from which training data will be extracted. 
#' Training data will be extracted from each layer in a \code{RasterBrick} or 
#' \code{RasterStack}.
#' @param polys a \code{SpatialPolygonsDataFrame} with training data polygons
#' @param classcol the name of the column containing the response variable (for 
#' example the land cover type of each pixel)
#' @return A data.frame
#' with the training data. Each row will contain the response (the column 
#' chosen by \code{classcol}) as the first column, with the remaining columns 
#' containing the values at that location of each band in the raster stack.  
#' @examples
#' L5TSR_1986 <- stack(system.file('extdata/L5TSR_1986.dat', package='teamr'))
#' data(L5TSR_1986_2001_training)
#' train_data <- extract_training_data(L5TSR_1986, L5TSR_1986_2001_training, 
#'                                     "t1_class")
extract_training_data <- function(x, polys, classcol, trainfrac=.7) {
    if (projection(x) != projection(polys)) {
        stop('Coordinate systems do not match')
    }
    # Convert classcol from the name of the column to an index
    classcol <- grep(paste0('^', classcol, '$'), names(polys))
    if (length(classcol) == 0) {
        stop(paste(classcol, 'not found in polys'))
    }
    pixels <- extract(x, polys, small=TRUE, df=TRUE)
    # Add column with class labels
    polys$ID <- seq(1, nrow(polys))
    pixels <- cbind(y=polys@data[match(pixels$ID, polys$ID), classcol], 
                    pixels)
    pixels <- pixels[!(names(pixels) == 'ID')]
    return(pixels)
}
