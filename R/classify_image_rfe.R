#' Runs an SVM image classification with recursive feature selection
#'
#' Uses the \code{rfe} function from the \code{caret} package to run a a 
#' recursive feature selection to select the best features to include in the
#' SVM model.
#'
#' @export
#' @import caret
#' @param x a \code{Raster*} image with the predictor layer(s) for the 
#' classification
#' @param train_data a data table with a column labeled 'y' with the observed 
#' classes, and one or more columns with the values of predictor(s) at each 
#' location.
#' @param classProbs whether to also calculate and return the probabilities of 
#' membership for each class
#' @param use_training_flag indicates whether to exclude data flagged as 
#' testing data when training the classifier. For this to work the input 
#' train_data \code{data.frame} must have a column named 'Training' that 
#' indicates, for each pixel, whether that pixel is a training pixel (coded as 
#' TRUE) or testing pixel (coded as FALSE).
#' @param tune_length the number of levels of each parameter that should be 
#' tried by \code{train} when training the classifier
#' @param train_control default is NULL (reasonable value will be set 
#' automatically).  For details see \code{\link{trainControl}}.
#' @param tune_grid the training grid to be used for training the classifier.  
#' Must be a \code{data.frame} with two columns: ".sigma" and ".C".
#' @param notify notifier to use (defaults to \code{print} function). See the 
#' \code{notifyR} package for one way of sending notifications from R. The 
#' \code{notify} function should accept a string as the only argument.
#' @return a list with 3 elements: the trained classifier, the predicted classes 
#' \code{RasterLayer} and the class probabilities \code{RasterBrick}
#' @details processing can be done in parallel using all available CPUs through 
#' the use of the cluster facilities in the \code{spatial.tools} package. To 
#' enable clustering, call \code{sfQuickInit} before running 
#' \code{classify_image}. To stop the cluster when finished, call 
#' \code{sfQuickStop}.
#'
#' @examples
#' \dontrun{
#' # Don't run long example
#' set.seed(0)
#' train_data_1986 <- extract_training_data(L5TSR_1986, 
#'                                          polys=L5TSR_1986_2001_training,
#'                                          class_col="class_1986", 
#'                                          training=.7)
#' classified_LT5SR_1986 <- classify_image(L5TSR_1986, train_data_1986)
#'
#' classified_LT5SR_1986$svmRFE
#' plot(classified_LT5SR_1986$pred_classes)
#' plot(classified_LT5SR_1986$pred_probs)
#' summary(accuracy(classified_LT5SR_1986$svmRFE))
#' }
classify_image_rfe <- function(x, train_data, classProbs=TRUE, 
                           use_training_flag=TRUE, tune_length=8, 
                           train_control=NULL, tune_grid=NULL, notify=print) {

    cl <- options('rasterClusterObject')[[1]]
    inparallel <- FALSE
    if (!is.null(cl)) {
        if (!require(doSNOW)) {
            warning('Cluster object found, but "doSNOW" package is required to run training in parallel. Running sequentially.')
        } else {
            registerDoSNOW(cl)
            inparallel <- TRUE
        }
    }

    notify('Training classifier...')
    if (is.null(train_control)) {
        train_control <- trainControl(method="repeatedcv",
                                      repeats=5, 
                                      classProbs=classProbs)
    }
    if (use_training_flag) {
        if (!('Training' %in% names(train_data))) {
            stop('when use_training_flag is TRUE, train_data must have a "Training" column')
        }
        train_data$x <- train_data$x[train_data$Training]
        train_data$y <- train_data$y[train_data$Training]
    }

    # This recursive feature elimination procedure follows Algorithm 19.5 in 
    # Kuhn and Johnson 2013
    svmFuncs <- caretFuncs
    
    # First center and scale
    normalization <- preProcess(train_data$x, method='range')
    scaled_predictors <- predict(normalization, train_data$x)
    scaled_predictors <- as.data.frame(scaled_predictors)
    subsets <- c(1:ncol(scaled_predictors))
    ctrl <- rfeControl(method="repeatedcv",
                       repeats=5,
                       verbose=TRUE,
                       functions=svmFuncs)

    svmRFE <- rfe(x=scaled_predictors,
                  y=train_data$y,
                  sizes=subsets,
                  metric="ROC",
                  rfeControl=ctrl,
                  method="svmRadial",
                  tuneLength=tune_length,
                  trControl=train_control)
                  
    notify('Predicting classes...')
    n_classes <- length(levels(svmRFE))
    if (inparallel) {
        pred_classes <- clusterR(x, predict, args=list(svmRFE))
    } else {
        pred_classes <- predict(x, svmRFE)
    }
    names(pred_classes) <- 'cover'

    if (classProbs) {
        notify('Calculating class probabilities...')
        if (inparallel) {
            pred_probs <- clusterR(x, predict, args=list(svmRFE, type="prob", 
                                                         index=c(1:n_classes)))
        } else {
            pred_probs <- predict(x, svmRFE, type="prob", index=c(1:n_classes))
        }
        return(list(model=svmRFE, pred_classes=pred_classes, 
                    pred_probs=pred_probs))
    } else {
        return(list(model=svmRFE, pred_classes=pred_classes))
    }
}
