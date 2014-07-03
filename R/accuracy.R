#' A class for representing accuracy assessment results
#' @exportClass accuracy
#' @rdname accuracy-class
#' @slot ct a simple sample contingency table
#' @slot pop_ct a population contingency table (if \code{pop} was provided - 
#' see \code{\link{accuracy}})
#' @slot Q quantity disagreement
#' @slot A allocation disagreement
#' @slot n_test the number of samples
#' @slot pop the population of each class as a numeric
#' @import methods
setClass('accuracy', slots=c(ct='table', pop_ct='table', Q='numeric', 
                             A='numeric', n_test='numeric',
                             pop='numeric')
)

#' @export
#' @method summary accuracy
summary.accuracy <- function(object, ...) {
    obj = list()
    obj[['class']] <- class(object)
    obj[['Q']] <- object@Q
    obj[['A']] <- object@A
    obj[['ct']] <- object@ct
    obj[['pop_ct']] <- object@pop_ct
    obj[['n_test']] <- object@n_test
    margined_pop_ct <- .add_ct_margins(object@pop_ct)
    obj[['overall_acc']] <- margined_pop_ct[length(margined_pop_ct)]
    class(obj) <- 'summary.accuracy'
    obj
}

#' @export
#' @method print summary.accuracy
print.summary.accuracy <- function(x, ...) {
    cat(paste('Object of class "', x[['class']], '"\n', sep = ''))
    cat('\n')
    cat(paste('Testing samples:\t', x[['n_test']], '\n', sep = ''))
    cat('\n')
    cat('Sample contingency table:\n')
    print(.add_ct_margins(x[['ct']]))
    cat('\n')
    cat('Population contingency table:\n')
    print(.add_ct_margins(x[['pop_ct']]))
    cat('\n')
    cat(paste("Overall accuracy:\t", round(x[['overall_acc']], digits=4), "\n", sep = ""))
    cat('\n')
    cat(paste('Quantity disagreement:\t\t', round(x[['Q']], digits=4), '\n', sep = ''))
    cat(paste('Allocation disagreement:\t', round(x[['A']], digits=4), '\n', sep = ''))
    invisible(x)
}

#' @export
#' @method print accuracy
print.accuracy <- function(x, ...) {
    print(summary(x, ...))
}

setMethod("show", signature(object="accuracy"), function(object) print(object))

#' A class for error adjusted class areas
#'
#' @seealso \code{\link{adj_areas}}.
#' @import methods
#' @exportClass error_adj_area
setClass('error_adj_area', slots=c(adj_area_mat='matrix'))

#' Calculated adjusted class areas for an image classification
#'
#' Calculates the adjusted areas of each class in an image after taking account 
#' of omission and commission errors. For unbiased adjustments, error rates 
#' should be calculated using a population sample matrix (see 
#' \code{\link{accuracy}}.
#'
#' Standard errors for the adjusted areas are calculated as in Olofsson et al.  
#' (2013).
#' @export adj_areas
#' @param x an \code{accuracy} object or a list of populations as a 
#' \code{numeric}
#' @param y missing, or a contingency table
#' @references Olofsson, P., G. M. Foody, S. V. Stehman, and C. E. Woodcock.  
#' 2013. Making better use of accuracy data in land change studies: Estimating 
#' accuracy and area and quantifying uncertainty using stratified estimation.  
#' Remote Sensing of Environment 129:122-131.
setGeneric("adj_areas", function(x, y) standardGeneric("adj_areas"))

#' @rdname adj_areas
#' @aliases adj_areas,numeric,table-method
setMethod("adj_areas", signature(x="numeric", y="table"),
function(x, y) {
    pop <- x
    ct <- y
    Wi <- pop / sum(pop)
    adj_area_est <- sum(pop) * colSums(Wi * (ct / rowSums(ct)))
    # Calculate standard errors of the proportions
    std_err_p <- sqrt(colSums(Wi^2 *
                              (((ct / rowSums(ct))*(1 - ct / rowSums(ct))) /
                               (rowSums(ct) - 1))))
    # Now calculate standard error of adjusted area estimate
    std_err_area <- sum(pop) * std_err_p
    adj_area_mat <- cbind(pop, adj_area_est, std_err_area, 1.96 * std_err_area)
    adj_area_mat <- round(adj_area_mat, 0)
    dimnames(adj_area_mat)[[1]] <- dimnames(ct)[[1]]
    dimnames(adj_area_mat)[[2]] <- c('Mapped area', 'Adj. area', 'S.E.', 
                                     '1.96 * S.E.')
    return(new('error_adj_area', adj_area_mat=adj_area_mat))
})

#' @rdname adj_areas
#' @aliases adj_areas,numeric,matrix-method
setMethod("adj_areas", signature(x="numeric", y="matrix"),
function(x, y) {
    class(y) <- "table"
    adj_areas(x, y)
})

#' @rdname adj_areas
#' @aliases adj_areas,numeric,missing-method
setMethod("adj_areas", signature(x="accuracy", y='missing'),
function(x) {
    pop <- x@pop
    ct <- x@ct
    adj_areas(pop, ct)
})

setMethod("show", signature(object="error_adj_area"),
function(object) {
    cat('Object of class: error_adj_area\n')
    cat('Accuracy-adjusted area table:\n')
    print(object@adj_area_mat)
})


plot.error_adj_area <- function(x, ...) {
    classes <- dimnames(x@adj_area_mat)[[1]]
    areas <- x@adj_area_mat[, 2]
    se <- x@adj_area_mat[, 3]
    plt_data <- data.frame(x=classes, y=areas, se=se)
    y <- NULL # Fix for R CMD check
    ggplot(plt_data, aes(x, y)) + geom_bar(stat="identity") + 
        geom_errorbar(aes(ymin=y - 1.96 * se, ymax=y + 1.96 * se), width=.25) +
        xlab("Class") + ylab("Area")
}

.calc_pop_ct <- function(ct, pop) {
    # Below uses the notation of Pontius and Millones (2011)
    nijsum <- matrix(rowSums(ct), nrow=nrow(ct), ncol=ncol(ct))
    Ni <- matrix(pop, nrow=nrow(ct), ncol=ncol(ct))
    # pop_ct is the population contigency table
    pop_ct <- (ct / nijsum) * (Ni / sum(pop))
    dimnames(pop_ct)[[1]] <- dimnames(ct)[[1]]
    dimnames(pop_ct)[[2]] <- dimnames(ct)[[2]]
    class(pop_ct) <- 'table'
    return(pop_ct)
}

.calc_Q <- function(pop_ct) {
    # Calculate quantity disagreement (Pontius and Millones, 2011, eqns 2-3)
    qg_mat = abs(rowSums(pop_ct) - colSums(pop_ct))
    return(sum(qg_mat) / 2)
}

.calc_A <- function(pop_ct) {
    # Calculate allocation disagreement (Pontius and Millones, 2011, eqns 4-5)
    diag_indices <- which(diag(nrow(pop_ct)) == TRUE)
    ag_mat = 2 * apply(cbind(rowSums(pop_ct) - pop_ct[diag_indices],
                             colSums(pop_ct) - pop_ct[diag_indices]), 1, min)
    return(sum(ag_mat) / 2)
}

# Adds margins to contingency table
.add_ct_margins <- function(ct, digits=4) {
    # For user's, producer's, and overall accuracy formulas, see Table 
    # 21.3 in Foody, G.M., Stehman, S.V., 2009. Accuracy Assessment, in: 
    # Warner, T.A., Nellis, M.D., Foody, G.M. (Eds.), The SAGE Handbook of 
    # Remote Sensing. SAGE.
    diag_indices <- which(diag(nrow(ct)) == TRUE)
    users_acc <- ct[diag_indices] / colSums(ct)
    prod_acc <- ct[diag_indices] / rowSums(ct)
    overall_acc <- sum(ct[diag_indices]) / sum(ct)
    ct <- addmargins(ct)
    dimnames(ct)[[1]][nrow(ct)] <- "Total"
    dimnames(ct)[[2]][nrow(ct)] <- "Total"
    ct <- rbind(ct, Producers=c(users_acc, NA))
    ct <- cbind(ct, Users=c(prod_acc, NA, overall_acc))
    ct <- round(ct, digits=digits)
    dimnames(ct) <- list(predicted=dimnames(ct)[[1]],
                         observed=dimnames(ct)[[2]])
    class(ct) <- 'table'
    return(ct)
}

#' Calculate statistics summarizing classification accuracy
#'
#' Calculates a contingency table and various statistics for use in image 
#' classification accuracy assessment and map comparison. Contingency table 
#' includes user's, producer's, and overall accuracies for an image 
#' classification, and quantity disagreement \code{Q} and allocation 
#' disagreement \code{A}. \code{Q} and \code{A} are calculated based on Pontius 
#' and Millones (2011). Standard errors for 95 percent confidence intervals for 
#' the user's, producer's and overall accuracies are calculated as in Foody and 
#' Stehman (2009) Table 21.3. To avoid bias due to the use of a sample 
#' contingency table, the contingency table will be converted to a population 
#' contingency table if the variable 'pop' is provided. For an accuracy 
#' assessment using testing data from a simple random sample, 'pop' does not 
#' need to be provided (see Details).
#'
#' \code{x} can be one of:
#' \enumerate{
#'
#'   \item A prediction model as output from one of the \code{teamlucc} 
#'   \code{classify} functions. If \code{x} is a model, and testing data 
#'   is included in the model, \code{pop} and \code{test_data} can both be 
#'   missing, and accuracy will still run (though the output will in this case 
#'   be biased unless the testing data is from a simple random sample). If 
#'   \code{x} is a \code{RasterLayer}, then \code{test_data} must be supplied.
#'
#'   \item A \code{RasterLayer} with a predicted map.
#' }
#'
#' \code{test_data} can be one of:
#' \enumerate{
#'   \item \code{NULL}. If test_data is \code{NULL}, \code{accuracy} will try to use 
#'         testing data included in \code{x}. This will only work if \code{x}
#'         is a model of class \code{train} from the \code{caret} package, and 
#'         if the model was run using the one of the \code{teamlucc} 
#'         \code{classify} functions.
#'
#'   \item A \code{SpatialPolygonsDataFrame} object, in which case \code{accuracy} 
#'         will extract the predicted classes within each polygon from \code{x}.  
#'         This will only work if \code{x} is a \code{RasterLayer}.
#'
#'   \item A \code{pixel_data} object, in which case \code{accuracy} will use the 
#'         included \code{training_flag} indicator to separate testing and 
#'         training data.
#' }
#'
#' \code{pop} can be one of:
#' \enumerate{
#'   \item NULL, in which case the sample frequencies will be used as estimates 
#'         of the population frequencies of each class.
#'
#'   \item A list of length equal to the number of classes in the map giving 
#'         the total number of pixels in the population for each class.
#'
#'   \item A predicted cover map from as a \code{RasterLayer}, from which the 
#'         class frequencies will be tabulated and used as the population 
#'         frequencies.
#' }
#' @export accuracy
#' @param x either a classification model with a \code{predict} method or a 
#' \code{RasterLayer} (see Details)
#' @param test_data a \code{link{pixel_data}} object, 
#' \code{SpatialPolygonsDataFrame}, or NULL (see Details).
#' @param pop A \code{RasterLayer}, \code{numeric} of length equal to the 
#' number of clasess, or NULL (see Details).
#' @param class_col required if \code{test_data} is a 
#' \code{SpatialPolygonsDataFrame}. Defines the name of the column containing 
#' the observed cover class IDs
#' @param reclass_mat a reclassification matrix to be used in the case of a 
#' model fit by \code{classify} with the \code{do_split} option selected
#' @return \code{\link{accuracy-class}} instance
#' @references Pontius, R. G., and M. Millones. 2011. Death to Kappa: birth of 
#' quantity disagreement and allocation disagreement for accuracy assessment.  
#' International Journal of Remote Sensing 32:4407-4429.
#'
#' Olofsson, P., G. M. Foody, S. V. Stehman, and C. E. Woodcock.  2013. Making 
#' better use of accuracy data in land change studies: Estimating accuracy and 
#' area and quantifying uncertainty using stratified estimation.  Remote 
#' Sensing of Environment 129:122-131.
#'
#' Foody, G.M., Stehman, S.V., 2009. Accuracy Assessment, in: Warner, T.A., 
#' Nellis, M.D., Foody, G.M. (Eds.), The SAGE Handbook of Remote Sensing. SAGE.
#' @examples
#' \dontrun{
#' train_data <- get_pixels(L5TSR_1986, L5TSR_1986_2001_training, "class_1986", 
#'                          training=.6)
#' model <- train_classifier(train_data)
#' accuracy(L5TSR_1986_rfmodel)
#' }
setGeneric("accuracy", function(x, test_data, pop, class_col, reclass_mat) 
           standardGeneric("accuracy"))

#' @rdname accuracy
#' @aliases accuracy,train,ANY,ANY,missing,ANY-method
setMethod("accuracy", signature(x="train", test_data="ANY", pop="ANY", class_col="missing", reclass_mat="ANY"),
    function(x, test_data, pop, class_col, reclass_mat) {
        if (missing(test_data)) {
            test_data <- x$trainingData
            names(test_data)[names(test_data) == '.outcome'] <- 'y'
        } else {
            test_data <- cbind(y=test_data@y, 
                               test_data@x,
                               training_flag=test_data@training_flag)
        }
        if (!('training_flag' %in% names(test_data))) {
            warning('no training_flag variable found - assuming none of "test_data" was used for model training')
        } else if (sum(test_data$training_flag == 1) == length(test_data$training_flag)) {
            stop('cannot conduct accuracy assessment without independent testing data')
        }
        test_data <- test_data[!test_data$training_flag, ]
        complete_rows <- complete.cases(test_data)
        if (sum(complete_rows) != nrow(test_data)) {
            warning(paste('ignored', nrow(test_data) - sum(complete_rows), 
                          'rows because of missing data'))
            test_data <- test_data[complete.cases(test_data), ]
        }
        predicted <- predict(x, test_data)
        observed <- test_data$y

        calc_accuracy(predicted, observed, pop, reclass_mat)
    }
)

#' @rdname accuracy
#' @aliases accuracy,RasterLayer,pixel_data,ANY,missing,ANY-method
setMethod("accuracy", signature(x="RasterLayer", test_data="pixel_data", pop="ANY", class_col="missing", reclass_mat="ANY"),
    function(x, test_data, pop, class_col, reclass_mat) {
        if (all(test_data@training_flag == 1)) {
            stop('cannot conduct accuracy assessment without independent testing data')
        } else if (all(test_data@training_flag == 0)) {
            # All the test_data is for testing
            predicted <- extract(x, test_data@polys, small=TRUE, df=TRUE)[, 2]
            observed <- test_data@y
        } else {
            # Mix of testing and validation data
            predicted <- extract(x, test_data@polys[!test_data@training_flag], 
                                 small=TRUE, df=TRUE)[, 2]
            observed <- test_data@y[!test_data@training_flag]
        }
        predicted <- factor(predicted, labels=levels(observed))
        calc_accuracy(predicted, observed, pop, reclass_mat)
    }
)

#' @rdname accuracy
#' @aliases accuracy,RasterLayer,SpatialPolygonsDataFrame,ANY,character,ANY-method
setMethod("accuracy", signature(x="RasterLayer", test_data="SpatialPolygonsDataFrame", pop="ANY", class_col="character", reclass_mat="ANY"),
    function(x, test_data, pop, class_col, reclass_mat) {
        ext <- get_pixels(x, test_data, class_col=class_col)
        # Since x is the predicted image, the output of get_pixels gives 
        # the predicted value in slot x, and the observed value in slot y.  
        # However x is converted to a numeric from a factor, so it needs to be 
        # converted back to a factor with the same levels as y.
        observed <- ext@y
        predicted <- factor(ext@x[, ], labels=levels(ext@y))
        calc_accuracy(predicted, observed, pop, reclass_mat)
    }
)

calc_accuracy <- function(predicted, observed, pop, reclass_mat) {
    if (!missing(reclass_mat)) {
        stop('reclass_mat not yet supported')
    }

    # ct is the sample contigency table
    ct <- table(predicted, observed)

    if (missing(pop)) {
        warning('pop was not provided - assuming sample frequencies equal population frequencies')
        pop <- rowSums(ct)
    } else if (class(pop) == 'RasterLayer') {
        pop <- freq(pop, useNA='no')[, 2]
        if (length(pop) != nrow(ct)) {
            stop('number of classes in pop must be equal to nrow(ct)')
        }
    } else if (class(pop) %in% c('integer', 'numeric')) {
        if (length(pop) != nrow(ct)) {
            stop('length(pop) must be equal to number of classes in the predicted data')
        }
    } else { 
        stop('pop must be a numeric vector or integer vector of length equal to the number of classes in x, or a RasterLayer, or NULL')
    }

    pop_ct <- .calc_pop_ct(ct, pop)
    Q <- .calc_Q(pop_ct)
    A <- .calc_A(pop_ct)

    return(new("accuracy", ct=ct, pop_ct=pop_ct, Q=Q, A=A, 
               n_test=length(observed), pop=pop))
}
