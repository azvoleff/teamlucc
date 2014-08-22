#' A class for representing training or testing data
#'
#' Used to represent training data for a machine learning classifier for image 
#' classificaion, or testing data used for testing a classification.
#'
#' @exportClass pixel_data
#' @rdname pixel_data-class
#' @aliases pixel_data
#' @slot x a \code{data.frame} of independent variables (usually pixel values)
#' @slot y a \code{data.frame} of the dependent variable (usually land cover 
#' classes)
#' @slot pixel_src a data.frame used to link pixels in \code{x} and \code{y} to 
#' an input polygon
#' @slot training_flag a binary vector of length equal to \code{nrow(x)} 
#' indicating each row in x should be used in training (TRUE) or in testing 
#' (FALSE)
#' @slot polys a \code{SpatialPolygonsDataFrame} of the polygons used to choose 
#' the pixels in \code{x} and \code{y}.
#' @import methods
#' @importFrom sp SpatialPolygonsDataFrame
setClass('pixel_data', slots=c(x='data.frame', y='factor', 
                               pixel_src='data.frame', training_flag='logical', 
                               polys='SpatialPolygonsDataFrame')
)

#' @export
#' @method summary pixel_data
summary.pixel_data <- function(object, ...) {
    obj = list()
    obj[['class']] <- class(object)
    obj[['n_classes']] <- nlevels(object)
    obj[['n_sources']] <- length(unique(object@polys$src))
    obj[['n_polys']] <- nrow(object@polys)
    obj[['n_pixels']] <- nrow(object@x)
    training_df <- data.frame(y=object@y,
                              pixel_src=src_name(object), 
                              training_flag=object@training_flag)
    y=pixel_src=training_flag=NULL # Keep R CMD CHECK happy
    class_stats <- summarize(group_by(training_df, y),
                             n_polys=length(unique(pixel_src)),
                             n_train_pixels=sum(training_flag),
                             n_test_pixels=sum(!training_flag),
                             train_frac=round(sum(training_flag) / 
                                              length(training_flag), 2))
    names(class_stats)[names(class_stats) == 'y'] <- 'class'
    obj[['class_stats']]  <- class_stats
    obj[['n_training']] <- sum(object@training_flag == TRUE)
    obj[['n_testing']] <- sum(object@training_flag == FALSE)
    obj[['training_frac']] <- sum(object@training_flag == TRUE) / length(object@training_flag)
    class(obj) <- 'summary.pixel_data'
    obj
}

#' @export
#' @method print summary.pixel_data
print.summary.pixel_data <- function(x, ...) {
    cat(paste('Object of class "', x[['class']], '"\n', sep = ''))
    cat('\n')
    cat(paste('Number of classes:\t', x[['n_classes']], '\n', sep=''))
    cat(paste('Number of polygons:\t', x[['n_polys']], '\n', sep=''))
    cat(paste('Number of pixels:\t', x[['n_pixels']], '\n', sep=''))
    cat(paste('Number of sources:\t', x[['n_sources']], '\n', sep=''))
    cat('\n')
    cat('Training data statistics:\n')
    print(x[['class_stats']])
    cat('\n')
    cat(paste('Number of training samples:\t', x[['n_training']], '\n', sep=''))
    cat(paste('Number of testing samples:\t', x[['n_testing']], '\n', sep=''))
    cat(paste('Training fraction:\t\t', round(x[['training_frac']], 2), '\n', sep=''))
    invisible(x)
}

#' @export
#' @method length pixel_data
length.pixel_data <- function(x) {
    return(length(x@y))
}

#' @export
#' @method levels pixel_data
levels.pixel_data <- function(x) {
    return(levels(x@y))
}

#' @export
#' @method print pixel_data
print.pixel_data <- function(x, ...) {
    print(summary(x, ...))
}

#' @export
#' @importFrom maptools spRbind
#' @method rbind pixel_data
rbind.pixel_data <- function(x, ...) {
    for (item in c(...)) {
        x@x <- rbind(x@x, item@x)
        x@y <- factor(c(as.character(x@y), as.character(item@y)))
        x@pixel_src <- rbind(x@pixel_src, item@pixel_src)
        x@training_flag <- c(x@training_flag, item@training_flag)
        if (any(row.names(x@polys) %in% row.names(item@polys)))
            stop('training polygon IDs are not unique - are src_names unique?')
        x@polys <- spRbind(x@polys, item@polys)
    }
    return(x)
}

#' Extract part of pixel_data class
#'
#' @method [ pixel_data
#' @rdname extract-methods
#' @param x a \code{pixel_data} object
#' @param i a class or list of classes to extract
#' @param j unused
#' @param ... additional arguments (none implemented)
setMethod("[", signature(x="pixel_data", i='character', j="ANY"),
function(x, i, j, ...) {
    if (!(i %in% levels(x@y))) {
        stop(paste0('"', i, '"', ' is not a class in this pixel_data object'))
    }
    sel_rows <- x@y %in% i
    used_polys <- which(paste(x@polys@data$src, x@polys@data$ID) %in% 
                        with(x@pixel_src[sel_rows, ], paste(src, ID)))
    initialize(x, x=x@x[sel_rows, ], y=x@y[sel_rows], 
               pixel_src=x@pixel_src[sel_rows, ], 
               training_flag=x@training_flag[sel_rows], 
               polys=x@polys[used_polys, ])
})

setMethod("show", signature(object="pixel_data"), function(object) 
          print(object))

#' Subsample a pixel_data object
#'
#' @export subsample
#' @param x a \code{pixel_data} object
#' @param size either 1) a number from 0 to 1, indicating \code{size} is the 
#' fraction of pixels to sample, or 2) a number greater than 1, in which case 
#' \code{size} is the number of pixels to sample. Size applies per strata, if 
#' stratification is chosen.
#' @param strata whether to draw samples from within individual classes, nested 
#' within source polygons (\code{strata='sources'}), or from within individual 
#' classes alone (\code{strata='classes'})
#' @param type whether to subsample training data (\code{type='training'}) or 
#' testing data (\code{type='testing'}). Whichever type is chosen, the other 
#' type will be left untouched (for example, if \code{type='testing'}, the 
#' training data will not be changed).
#' @param flag whether to swap training flag on sampled data (for example, flag 
#' sampled training data as testing data, if \code{flag=TRUE} and 
#' \code{type='training'}) or remove sampled data from dataset entirely 
#' (\code{flag=FALSE}).
#' @param classes specifies which classes to sample, defaults to all classes in 
#' \code{x}
#' @rdname subsample
#' @aliases subsample,pixel_data-method
setGeneric("subsample", function(x, size, strata="sources", type="training", 
                                 flag=TRUE, classes=levels(x@y))
    standardGeneric("subsample")
)

#' @rdname subsample
#' @aliases subsample,pixel_data,numeric-method
#' @importFrom dplyr group_by sample_frac
setMethod("subsample", signature(x="pixel_data", size="numeric"),
function(x, size, strata, type, flag, classes) {
    row_IDs <- data.frame(y=x@y,
                          pixel_src=paste(x@pixel_src$src, x@pixel_src$ID),
                          row_num=seq(1, length(x@y)))
    stopifnot(size > 0)
    stopifnot(strata %in% c("sources", "classes"))
    stopifnot(type %in% c("training", "testing"))
    if (type == "training") {
        row_IDs <- row_IDs[x@training_flag, ]
    } else {
        row_IDs <- row_IDs[!x@training_flag, ]
    }
    row_IDs <- row_IDs[row_IDs$y %in% classes, ]
    stopifnot(nrow(row_IDs) > 1)
    if (strata == "sources") {
        y=pixel_src=NULL
        row_IDs <- group_by(row_IDs, y, pixel_src)
    } else if (strata == "classes") {
        row_IDs <- group_by(row_IDs, y)
    }
    if (size < 1) {
        samp_rows <- dplyr:::sample_frac.grouped_df(row_IDs, size)$row_num
    } else {
        samp_rows <- dplyr:::sample_n.grouped_df(row_IDs, size)$row_num
    }
    if (flag) {
        if (type == 'testing') {
            x@training_flag[samp_rows] <- TRUE
        } else if (type == 'training') {
            x@training_flag[samp_rows] <- FALSE
        }
    } else {
        x@x <- x@x[samp_rows, ]
        x@y <- x@y[samp_rows]
        x@training_flag <- x@training_flag[samp_rows]
        x@pixel_src <- x@pixel_src[samp_rows, ]
    }
    return(x)
})

#' Get or set training_flag for a pixel_data object
#'
#' @export training_flag
#' @param x a \code{pixel_data} object
#' @param classes specifies a subset of classes in \code{x}
#' @aliases training_flag,pixel_data-method
setGeneric("training_flag", function(x, classes=levels(x@y)) {
    standardGeneric("training_flag")
})

#' @rdname training_flag
setMethod("training_flag", signature(x="pixel_data"),
function(x, classes) {
    if (identical(classes, levels(x@y))) {
        return(x@training_flag)
    } else {
        return(x@training_flag[x@y %in% classes])
    }
})

#' @rdname training_flag
#' @export
#' @param value training flag to assign for pixels in \code{x}
setGeneric("training_flag<-", function(x, classes=levels(x@y), value) {
    standardGeneric("training_flag<-")
})

#' @rdname training_flag
setMethod("training_flag<-", signature(x="pixel_data"),
function(x, classes=levels(x@y), value) {
    if (identical(classes, levels(x@y))) {
        # More efficiently handle special case of reassigning flags for all 
        # classes in x.
        if (length(value) == 1) value <- rep(value, length(x@training_flag))
        stopifnot(length(value) == length(x@training_flag))
        x@training_flag <- value
        return(x)
    } else {
        sel_rows <- which(x@y %in% classes)
        if (length(value) == 1) value <- rep(value, length(sel_rows))
        stopifnot(length(value) == length(sel_rows))
        x@training_flag[sel_rows] <- value
        return(x)
    }
})

#' Get number of testing pixels in a pixel_data object
#'
#' @export n_test
#' @param x a \code{pixel_data} object
#' @param classes specifies a subset of classes in \code{x}
#' @aliases n_test,pixel_data-method
setGeneric("n_test", function(x, classes=levels(x@y)) {
    standardGeneric("n_test")
})

#' @rdname n_test
setMethod("n_test", signature(x="pixel_data"),
function(x, classes) {
    if (identical(classes, levels(x@y))) {
        return(sum(!x@training_flag))
    } else {
        return(sum(!x@training_flag[x@y %in% classes]))
    }
})

#' Get number of training pixels in a pixel_data object
#'
#' @export n_train
#' @param x a \code{pixel_data} object
#' @param classes specifies a subset of classes in \code{x}
#' @aliases n_train,pixel_data-method
setGeneric("n_train", function(x, classes=levels(x@y)) {
    standardGeneric("n_train")
})

#' @rdname n_train
setMethod("n_train", signature(x="pixel_data"),
function(x, classes) {
    if (identical(classes, levels(x@y))) {
        return(sum(x@training_flag))
    } else {
        return(sum(x@training_flag[x@y %in% classes]))
    }
})

#' Get or set src_name for a pixel_data object
#'
#' @export src_name
#' @param x a \code{pixel_data} object
#' @param classes specifies a subset of classes in \code{x}
#' @aliases src_name,pixel_data-method
setGeneric("src_name", function(x, classes=levels(x@y)) {
    standardGeneric("src_name")
})

#' @method src_name pixel_data
setMethod("src_name", signature(x="pixel_data"),
function(x, classes) {
    if (identical(classes, levels(x@y))) {
        return(paste0(x@pixel_src$src, '_', x@pixel_src$ID))
    } else {
        return(with(x@pixel_src[x@y %in% classes, ], paste0(src, '_', ID)))
    }
})

#' @export
#' @rdname src_name
#' @param value a new \code{src_name} to assign for pixels in \code{x}
setGeneric("src_name<-", function(x, value) standardGeneric("src_name<-"))

#' @rdname src_name
setMethod("src_name<-", signature(x="pixel_data"),
function(x, value) {
    if (length(value) == 1) {
        value <- rep(value, nrow(x@polys))
    } else if (length(value) != nrow(x@polys)) {
        stop('src_name must be equal to 1 or number of polygons in x')
    }
    old_full_polyID <- paste(x@polys$src, x@polys$ID)
    x@polys$src <- value
    row.names(x@polys) <- paste0(x@polys$src, '_', x@polys$ID)

    new_full_polyID <- paste(x@polys$src, x@polys$ID)

    poly_pixel_match <- match(paste(x@pixel_src$src, x@pixel_src$ID), 
                              old_full_polyID)

    x@pixel_src$src <- x@polys$src[poly_pixel_match]
    x@pixel_src$ID <- x@polys$ID[poly_pixel_match]

    return(x)
})

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
#' used in training (rows equal to 1) or in testing (rows equal to FALSE), or 
#' 2) a logical vector of length equal to length(polys), or 3) a number between 
#' 0 and 1 indicating the fraction of the polygons to be randomly selected for 
#'   use in training.
#' @param src name of this data source. Useful when gathering training 
#' data from multiple images.
#' @return a \code{link{pixel_data}} object
#' will contain the the @examples
#' set.seed(1)
#' train_data <- get_pixels(L5TSR_1986, L5TSR_1986_2001_training, "class_1986", 
#'                          training=.6)
get_pixels <- function(x, polys, class_col, training=1, src='none') {
    if (projection(x) != projection(polys)) {
        stop('Coordinate systems do not match')
    }
    stopifnot(length(src) == 1)
    if (tolower(class_col) == 'id') {
        stop('class_col cannot be named "ID" (case insensitive)')
    }
    # Convert class_col from the name of the column to an index
    class_colnum <- grep(paste0('^', class_col, '$'), names(polys))
    if (length(class_colnum) == 0) {
        stop(paste0('"', class_col, '" not found in polys'))
    }
    # This is displayed in the dataframe, and should never change
    polys$ID <- row.names(polys)
    polys$src <- src
    row.names(polys) <- paste0(polys$src, '_', polys$ID)
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

    poly_rows <- pixels$ID
    pixels <- pixels[!(names(pixels) == 'ID')]

    # Convert y classes to valid R variable names - if they are not valid R 
    # variable names, the classification algorithm may throw an error
    y <- factor(make.names(polys@data[poly_rows, class_colnum]))

    pixel_src <- data.frame(src=polys@data[poly_rows, ]$src,
                            ID=polys@data[poly_rows, ]$ID, 
                            stringsAsFactors=FALSE)

    return(new("pixel_data", x=pixels, y=y, pixel_src=pixel_src,
               training_flag=polys@data[poly_rows, ]$training_flag,
               polys=polys))
}
