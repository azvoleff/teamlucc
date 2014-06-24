#' Trains a classifier
#'
#' This function trains a Support Vector Machine (SVM) or Random Forest (RF) 
#' classifier for use in an image classification.
#'
#' For \code{type='svm'}, \code{tunegrid} must be a \code{data.frame} with two 
#' columns: ".sigma" and ".C". For \code{type='rf'}, must be a 
#' \code{data.frame} with one column: '.mtry'.
#'
#' @export
#' @import caret randomForest
#' @param train_data a \code{link{pixel_data}} object
#' @param type either "svm" (to fit a support vector machine) or "rf" (to fit a
#' random forest).
#' @param use_training_flag indicates whether to exclude data flagged as 
#' testing data when training the classifier. For this to work the input 
#' train_data \code{data.frame} must have a column named 'training_flag' that 
#' indicates, for each pixel, whether that pixel is a training pixel (coded as 
#' TRUE) or testing pixel (coded as FALSE).
#' @param train_control default is NULL (reasonable values will be set 
#' automatically).  For details see \code{\link{trainControl}}.
#' @param tune_grid the training grid to be used for training the classifier.  
#' See Details.
#' @param use_rfe whether to use Recursive Feature Extraction (RFE) as 
#' implemented in the \code{caret} package to select a subset of the input 
#' features to be used in the classification. NOT YET SUPPORTED.
#' @param notify notifier to use (defaults to \code{print} function). See the 
#' \code{notifyR} package for one way of sending notifications from R. The 
#' \code{notify} function should accept a string as the only argument.
#' @return a trained random forest model (as a \code{train} object from the  
#' \code{caret} package)
train_classifier <- function(train_data, type='rf', use_training_flag=TRUE, 
                             train_control=NULL, tune_grid=NULL,
                             use_rfe=FALSE) {
    stopifnot(type %in% c('svm', 'rf'))

    if (is.null(train_control)) {
        train_control <- trainControl(method="oob")
    }

    # Build the formula, excluding the training flag column (if it exists) from 
    # the model formula
    predictor_names <- names(train_data@x)
    model_formula <- formula(paste('y ~', paste(predictor_names,
                                                collapse=' + ')))

    if (use_rfe) {
        stop('recursive feature extraction not yet supported')
        notify('Training classifier using RFE...')
        # This recursive feature elimination procedure follows Algorithm 19.5 
        # in Kuhn and Johnson 2013
        svmFuncs <- caretFuncs
        # First center and scale
        normalization <- preProcess(train_data@x, method='range')
        scaled_predictors <- predict(normalization, train_data@x)
        scaled_predictors <- as.data.frame(scaled_predictors)
        subsets <- c(1:length(predictor_names))
        ctrl <- rfeControl(method="repeatedcv",
                           repeats=5,
                           verbose=TRUE,
                           functions=svmFuncs)
        # For the rfe modeling, extract the training data from the main
        # train_data dataset - no need to pass the testing data to rfe
        rfe_x <- scaled_predictors[train_data@training_flag, ]
        rfe_y <- train_data@y[train_data@training_flag, ]
        rfe_res <- rfe(x=rfe_x, rfe_y,
                       sizes=subsets,
                       metric="ROC",
                       rfeControl=ctrl,
                       method="svmRadial",
                       tuneLength=tune_length,
                       trControl=train_control,
                       tuneGrid=tune_grid)
        #TODO: Extract best model from rfe_res
    } else {
        rfe_res <- NULL
    }

    train_data <- cbind(y=train_data@y,
                        train_data@x,
                        training_flag=train_data@training_flag,
                        poly_src=train_data@pixel_src$src,
                        poly_ID=train_data@pixel_src$ID)

    if (type == 'rf') {
        model <- train(model_formula, data=train_data, method="rf",
                       preProc=c('range'), subset=train_data$training_flag,
                       trControl=train_control, tuneLength=10, 
                       tuneGrid=tune_grid, ntree=2000)
    } else if (type == 'svm') {
        model <- train(model_formula, data=train_data, method="svmRadial",
                       preProc=c('range'), subset=train_data$training_flag,
                       trControl=train_control, tuneLength=10, 
                       tuneGrid=tune_grid)
    } else {
        # should never get here
        stop("model type not recognized")
    }

    return(model)
}
