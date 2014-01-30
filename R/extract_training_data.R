#' Extract training data for use in a classification
#'
#' @export
#' @param x a \code{Raster*} object from which training data will be extracted.  
#' Training data will be extracted from each layer in a \code{RasterBrick} or 
#' \code{RasterStack}.
#' @param polys a \code{SpatialPolygonsDataFrame} with training data polygons
#' @param class_col the name of the column containing the response variable (for 
#' example the land cover type of each pixel)
#' @param training indicator of which polygons to use in training. Can be: 1) a 
#' string giving the name of a column indicating whether each polygon is to be 
#' used in training (column equal to TRUE) or in testing (column equal to 
#' FALSE), or 2) a logical vector of length equal to length(polys), or 3) a 
#' number between 0 and 1 indicating the fraction of the polygons to be 
#' randomly selected for use in training.
#' @return A data.frame
#' with the training data. Each row will contain the response (the column 
#' chosen by \code{class_col}) as the first column, with the remaining columns 
#' containing the values at that location of each band in the raster stack.  
#' @examples
#' set.seed(1)
#' train_data <- extract_training_data(L5TSR_1986, L5TSR_1986_2001_training, 
#'                                     "class_1986", training=.6)
extract_training_data <- function(x, polys, class_col, training=1) {
    if (projection(x) != projection(polys)) {
        stop('Coordinate systems do not match')
    }
    # Convert class_col from the name of the column to an index
    class_colnum <- grep(paste0('^', class_col, '$'), names(polys))
    if (length(class_colnum) == 0) {
        stop(paste0('"', class_col, '" not found in polys'))
    }
    polys$FID <- row.names(polys)
    if (is.character(training)) {
        # Handle case of having column name suppled as 'training'
        training_col_index <- grep(training, names(polys))
        if (length(training_col_index) == 0) {
            stop(paste0('"', training,  '" column not found in x'))
        } else if (!is.logical(polys[training_col_index])) {
            stop(paste0('"', training,  '" column must be a logical vector'))
        }
        training <- polys[, grep(training, names(polys))]
    } else if (is.numeric(training) && (length(training) == 1) &&
               (training >= 0) && (training <= 1)) {
        # Handle case of having fraction supplied as 'training'
        if ('Training' %in% names(polys)) {
            stop('"Training" column already present in polys')
        }
        if (training == 0) {
            # Handle training=0 separately to enable use of quantile function 
            # below.
            polys$Training <- FALSE
        } else {
            sample_strata <- function(x) {
                      rand_vals <- runif(length(x))
                      rand_vals <= quantile(rand_vals, training)
            }
            polys$Training <- unlist(tapply(polys@data$FID, 
                                            polys@data[class_colnum], 
                                            sample_strata))
        }
    } else if ((length(training) == length(polys)) && is.logical(training)) {
        # Handle case of having vector supplied as 'training'
        polys$Training <- training
    } else {
        stop('"training" must be a column name, vector of same length as polys, or length 1 numeric')
    }
    pixels <- extract(x, polys, small=TRUE, df=TRUE)
    poly_pixel_match <- match(pixels$ID, seq(1, nrow(polys)))
    pixels <- pixels[!(names(pixels) == 'ID')]

    # Convert y classes to valid R variable names - if they are not valid R 
    # variable names, the classification algorithm may throw an error
    y <- factor(make.names(polys@data[poly_pixel_match, class_colnum]))

    return(list(y=y,
                x=pixels,
                Poly_FID=polys@data[poly_pixel_match, ]$FID,
                Training=polys@data[poly_pixel_match, ]$Training))
}
