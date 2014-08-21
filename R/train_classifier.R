#' Train a random forest or SVM classifier
#'
#' This function trains a Support Vector Machine (SVM) or Random Forest (RF) 
#' classifier for use in an image classification.
#'
#' For \code{type='svm'}, \code{tunegrid} must be a \code{data.frame} with two 
#' columns: ".sigma" and ".C". For \code{type='rf'}, must be a 
#' \code{data.frame} with one column: '.mtry'.
#'
#' This function will run in parallel if a parallel backend is registered with 
#' \code{\link{foreach}}.
#'
#' @export
#' @import caret randomForest e1071 kernlab
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
#' @param factors a list of character vector giving the names of predictors 
#' (layer names from the images used to build \code{train_data}) that should be 
#' treated as factors, and specifying the levels of each factor. For example, 
#' \code{factors=list(year=c(1990, 1995, 2000, 2005, 2010))}.
#' @param ... additional arguments (such as \code{ntree} for random forest 
#' classifier) to pass to \code{train}
#' @return a trained model (as a \code{train} object from the \code{caret} 
#' package)
#' @examples
#' train_data <- get_pixels(L5TSR_1986, L5TSR_1986_2001_training, "class_1986", 
#'                          training=.6)
#' model <- train_classifier(train_data)
train_classifier <- function(train_data, type='rf', use_training_flag=TRUE, 
                             train_control=NULL, tune_grid=NULL,
                             use_rfe=FALSE, factors=list(), ...) {
    stopifnot(type %in% c('svm', 'rf'))

    predictor_names <- names(train_data@x)

    # Convert predictors in training data to factors as necessary
    stopifnot(length(factors) == 0 || all(names(factors) %in% predictor_names))
    stopifnot(length(unique(names(factors))) == length(factors))
    for (factor_var in names(factors)) {
        pred_index <- which(predictor_names == factor_var)
        train_data@x[, pred_index] <- factor(train_data@x[, pred_index], 
                                             levels=factors[[factor_var]])
    }

    # Build the formula, excluding the training flag column (if it exists) from 
    # the model formula
    model_formula <- formula(paste('y ~', paste(predictor_names, collapse=' + ')))

    if (use_rfe) {
        stop('recursive feature extraction not yet supported')
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
        if (is.null(train_control)) {
            train_control <- trainControl(method="oob", classProbs=TRUE)
        }
        model <- train(model_formula, data=train_data, method="rf",
                       subset=train_data$training_flag,
                       trControl=train_control, tuneGrid=tune_grid, ...)
    } else if (type == 'svm') {
        if (is.null(train_control)) {
            train_control <- trainControl(method="cv", classProbs=TRUE)
        }
        model <- train(model_formula, data=train_data, method="svmRadial",
                       preProc=c('center', 'scale'), subset=train_data$training_flag,
                       trControl=train_control, tuneGrid=tune_grid, ...)
    } else {
        # should never get here
        stop("model type not recognized")
    }

    return(model)
}
