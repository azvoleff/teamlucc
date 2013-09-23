#' Extract training data for use in a classification
#'
#' @export
#' @param x A \code{Raster*} object from which training data will be 
#' extracted. Training data will be extracted from each layer in a 
#' \code{RasterBrick} or \code{RasterStack}.
#' @param y A shapefile with class polygons.
#' @return A data.frame with the training data. Each row will contain the 
#' response ('y') as the first column, with the remaining columns containing 
#' the values at that location of each band in the raster stack. 
#' @examples
#' # TODO: Add example
#' y <- 'H:/Data/TEAM/VB/Vectors/VB_training_1986_037_LT5.shp'
#' img_file <- 'H:/Data/TEAM/VB/Rasters/Landsat/1986_037_LT5/proc/lndsr.LT50150531986037XXX18.bsq'
#' x <- stack(img_file, bands=as.integer(c(1, 2, 3, 4)))
#' extract_training_data(x, y)
extract_training_data <- function(x, y) {
    layer_name <- gsub('.shp$', '', basename(y), ignore.case=TRUE)
    class_polys <- readOGR(dirname(y), layer_name, stringsAsFactors=FALSE)
    if (projection(x) != projection(class_polys)) {
        stop('Coordinate systems do not match')
    }
    pixel_sets <- extract(x, class_polys, small=TRUE)
    # Label pixel sets
    train_data <- mapply(function(class_type, pixel_set) 
                         cbind("y"=rep(class_type, nrow(pixel_set)), 
                               pixel_set), class_polys$Poly_Type, pixel_sets)
    # Convert training data to a data.frame
    train_data <- data.frame(abind(train_data, along=1), row.names=NULL)
    return(train_data)
}
