#' Runs an image classification
#'
#' Currently only supports classification using a random forest
#'
#' @export
#' @import caret randomForest
#' @param x a \code{Raster*} image with the predictor layer(s) for the 
#' classification
#' @param train_data a \code{link{Training_data}} object
#' @param class_probs whether to also calculate and return the probabilities of 
#' membership for each class
#' @param use_training_flag indicates whether to exclude data flagged as 
#' testing data when training the classifier. For this to work the input 
#' train_data \code{data.frame} must have a column named 'training_flag' that 
#' indicates, for each pixel, whether that pixel is a training pixel (coded as 
#' TRUE) or testing pixel (coded as FALSE).
#' @param train_control default is NULL (reasonable value will be set 
#' automatically).  For details see \code{\link{trainControl}}.
#' @param tune_grid the training grid to be used for training the classifier.  
#' Must be a \code{data.frame} with two columns: ".sigma" and ".C".
#' @param use_rfe whether to use Recursive Feature Extraction (RFE) as 
#' implemented in the \code{caret} package to select a subset of the input 
#' features to be used in the classification
#' @param notify notifier to use (defaults to \code{print} function). See the 
#' \code{notifyR} package for one way of sending notifications from R. The 
#' \code{notify} function should accept a string as the only argument.
#' @return a list with 3 elements: the trained classifier, the predicted classes 
#' \code{RasterLayer} and the class probabilities \code{RasterBrick}
#' @details Processing can be done in parallel using all using the cluster 
#' facilities in the \code{spatial.tools} package. To enable clustering, call 
#' \code{beginCluster} before running \code{classify_image_rf}.  To stop the 
#' cluster when finished, call \code{endCluster}.
#' @examples
#' \dontrun{
#' # Don't run long example
#' set.seed(0)
#' train_data_1986 <- extract_observed(L5TSR_1986, 
#'                                          polys=L5TSR_1986_2001_training,
#'                                          class_col="class_1986", 
#'                                          training=.7)
#' classified_LT5SR_1986 <- classify_image_rf(L5TSR_1986, train_data_1986)
#'
#' classified_LT5SR_1986$model
#' plot(classified_LT5SR_1986$pred_classes)
#' plot(classified_LT5SR_1986$pred_probs)
#' accuracy(classified_LT5SR_1986$model)
#' }
classify_image_rf <- function(x, train_data, class_probs=TRUE, 
                           use_training_flag=TRUE, 
                           train_control=NULL, tune_grid=NULL, 
                           use_rfe=FALSE, notify=print) {
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

    if (is.null(train_control)) {
        train_control <- trainControl(classProbs=class_probs)
    }
    # Build the formula, excluding the training flag column (if it exists) from 
    # the model formula
    predictor_names <- names(train_data@x)
    model_formula <- formula(paste('y ~',
                                   paste(predictor_names, collapse=' + ')))

    if (use_rfe) {
        stop('rfe not yet supported')
    }
    # if (use_rfe) {
    #     notify('Training classifier using RFE...')
    #     # This recursive feature elimination procedure follows Algorithm 19.5 
    #     # in Kuhn and Johnson 2013
    #     svmFuncs <- caretFuncs
    #     # First center and scale
    #     normalization <- preProcess(train_data@x, method='range')
    #     scaled_predictors <- predict(normalization, train_data@x)
    #     scaled_predictors <- as.data.frame(scaled_predictors)
    #     subsets <- c(1:length(predictor_names))
    #     ctrl <- rfeControl(method="repeatedcv",
    #                        repeats=5,
    #                        verbose=TRUE,
    #                        functions=svmFuncs)
    #     # For the rfe modeling, extract the training data from the main # 
    #     train_data dataset - no need to pass the testing data to rfe
    #     rfe_x <- scaled_predictors[train_data@training_flag, ]
    #     rfe_y <- train_data@y[train_data@training_flag, ]
    #     rfe_res <- rfe(x=rfe_x, rfe_y,
    #                    sizes=subsets,
    #                    metric="ROC",
    #                    rfeControl=ctrl,
    #                    method="svmRadial",
    #                    tuneLength=tune_length,
    #                    trControl=train_control,
    #                    tuneGrid=tune_grid)
    #     #TODO: Extract best model from rfe_res
    # } else {
    #     rfe_res <- NULL
    # }

    train_data <- cbind(y=train_data@y,
                        train_data@x,
                        training_flag=train_data@training_flag,
                        poly_ID=train_data@poly_ID)

    notify('Training classifier...')
    model <- train(model_formula, data=train_data, method="rf",
                   preProc=c('range'), subset=train_data$training_flag,
                   trControl=train_control, tuneGrid=tune_grid)

    notify('Predicting classes...')
    n_classes <- length(levels(model))
    if (inparallel) {
        pred_classes <- clusterR(x, predict, args=list(model))
    } else {
        pred_classes <- predict(x, model)
    }
    names(pred_classes) <- 'cover'


    if (class_probs) {
        notify('Calculating class probabilities...')
        if (inparallel) {
            pred_probs <- clusterR(x, predict, args=list(model, type="prob", index=c(1:n_classes)))
        } else {
            pred_probs <- predict(x, model, type="prob", index=c(1:n_classes))
        }
    } else {
        pred_probs <- NULL
    }

    return(list(model=model, pred_classes=pred_classes, pred_probs=pred_probs))
}
