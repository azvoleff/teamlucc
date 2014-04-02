#' Class to represent training data for training machine learning algorithms 
#'
#' @import methods
#' @importFrom sp SpatialPolygonsDataFrame
#' @export
#' @name Training_data-class
setClass('Training_data', slots=c(x='data.frame', y='factor', 
                                  poly_ID='character', training_flag='logical', 
                                  polys='SpatialPolygonsDataFrame')
)

#' @importFrom plyr ddply summarize .
#' @S3method summary Training_data
summary.Training_data <- function(object, ...) {
    obj = list()
    obj[['class']] <- class(object)
    obj[['n_classes']] <- nlevels(object)
    obj[['n_polys']] <- length(unique(object@poly_ID))
    obj[['n_pixels']] <- nrow(object@x)
    training_df <- data.frame(y=object@y, poly_ID=object@poly_ID, 
                              training_flag=object@training_flag)
    y=poly_ID=training_flag=NULL # Keep R CMD CHECK happy
    class_stats <- ddply(training_df, .(class=y), summarize,
                         n_pixels=length(y),
                         n_polys=length(unique(poly_ID)),
                         train_frac=round(sum(training_flag) / length(training_flag), 2))
    obj[['class_stats']]  <- class_stats
    obj[['training_frac']] <- sum(object@training_flag==TRUE) / length(object@training_flag)
    class(obj) <- 'summary.Training_data'
    obj
}

#' @S3method print summary.Training_data
print.summary.Training_data <- function(x, ...) {
    cat(paste('Object of class "', x[['class']], '"\n', sep = ''))
    cat('\n')
    cat(paste('Number of classes:\t', x[['n_classes']], '\n', sep=''))
    cat(paste('Number of polygons:\t', x[['n_polys']], '\n', sep=''))
    cat(paste('Number of pixels:\t', x[['n_pixels']], '\n', sep=''))
    cat('\n')
    cat('Training data statistics:\n')
    print(x[['class_stats']])
    cat('\n')
    cat(paste('Training fraction:\t', round(x[['training_frac']], 2), '\n', sep=''))
    invisible(x)
}

#' @S3method levels Training_data
levels.Training_data <- function(x) {
    return(levels(x@y))
}

#' @S3method print Training_data
print.Training_data <- function(x, ...) {
    print(summary(x, ...))
}

#' Show a Training_data object
#'
#' @export
setMethod("show", signature(object="Training_data"), function(object) 
          print(object))

#' Extract observed data for use in a classification (training or testing)
#'
#' @export
#' @param x a \code{Raster*} object from which observed data will be extracted.  
#' The data will be extracted from each layer in a \code{RasterBrick} or 
#' \code{RasterStack}.
#' @param polys a \code{SpatialPolygonsDataFrame} with polygons, each of which 
#' has been assigned to a particular class (using the \code{class_col}
#' @param class_col the name of the column containing the response variable 
#' (for example the land cover type of each pixel)
#' @param training indicator of which polygons to use in training. Can be: 1) a 
#' string giving the name of a column indicating whether each polygon is to be 
#' used in training (column equal to TRUE) or in testing (column equal to 
#' FALSE), or 2) a logical vector of length equal to length(polys), or 3) a 
#' number between 0 and 1 indicating the fraction of the polygons to be 
#' randomly selected for use in training.
#' @return data.frame with the training data. Each row will contain the 
#' response (the column chosen by \code{class_col}) as the first column, with 
#' the remaining columns containing the values at that location of each band in 
#' the raster stack.
#' @examples
#' set.seed(1)
#' train_data <- extract_observed(L5TSR_1986, L5TSR_1986_2001_training, 
#'                                "class_1986", training=.6)
extract_observed <- function(x, polys, class_col, training=1) {
    if (projection(x) != projection(polys)) {
        stop('Coordinate systems do not match')
    }
    # Convert class_col from the name of the column to an index
    class_colnum <- grep(paste0('^', class_col, '$'), names(polys))
    if (length(class_colnum) == 0) {
        stop(paste0('"', class_col, '" not found in polys'))
    }
    polys$ID <- row.names(polys)
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
        if ('training_flag' %in% names(polys)) {
            stop('"training_flag" column already present in polys')
        }
        if (training == 0) {
            # Handle training=0 separately to enable use of quantile function 
            # below.
            polys$training_flag <- FALSE
        } else {
            sample_strata <- function(x) {
                rand_vals <- runif(length(x))
                rand_vals <= quantile(rand_vals, training)
            }
            polys$training_flag <- unlist(tapply(polys@data$ID, 
                                                 polys@data[class_colnum], 
                                                 sample_strata))
        }
    } else if ((length(training) == length(polys)) && is.logical(training)) {
        # Handle case of having vector supplied as 'training'
        polys$training_flag <- training
    } else {
        stop('"training" must be a column name, vector of same length as polys, or length 1 numeric')
    }
    pixels <- extract(x, polys, small=TRUE, df=TRUE)
    poly_pixel_match <- match(pixels$ID, seq(1, nrow(polys)))
    pixels <- pixels[!(names(pixels) == 'ID')]

    # Convert y classes to valid R variable names - if they are not valid R 
    # variable names, the classification algorithm may throw an error
    y <- factor(make.names(polys@data[poly_pixel_match, class_colnum]))

    return(new("Training_data", x=pixels, y=y, 
               poly_ID=polys@data[poly_pixel_match, ]$ID,
               training_flag=polys@data[poly_pixel_match, ]$training_flag,
               polys=polys))
}
