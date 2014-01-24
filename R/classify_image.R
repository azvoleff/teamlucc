#' Runs an image classification
#'
#' Currently only supports classification using a support vector machine (SVM).
#'
#' @export
#' @import caret
#' @param x a \code{Raster*} image with the predictor layer(s) for the 
#' classification
#' @param train_data a data table with a column labeled 'y' with the observed 
#' classes, and one or more columns with the values of predictor(s) at each 
#' location.
#' @param pred_classes_filename the filename to use to save the predicted 
#' classes \code{RasterLayer}. If 'NULL' (the default) a temporary file will be 
#' used if necessary (if the size of the output raster exceeds available 
#' memory).
#' @param pred_probs_filename the filename to use to save the class 
#' probabilities \code{RasterLayer} or \code{RasterBrick}. If 'NULL' (the 
#' default) a temporary file will be used if necessary (if the size of the 
#' output raster exceeds available memory).
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
#' @param train_grid the training grid to be used for training the classifier.  
#' Must be a \code{data.frame} with two columns: ".sigma" and ".C".
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
#'                                          classcol="class_1986", training=.7)
#' classified_LT5SR_1986 <- classify_image(L5TSR_1986, train_data_1986)
#'
#' classified_LT5SR_1986$model
#' plot(classified_LT5SR_1986$pred_classes)
#' plot(classified_LT5SR_1986$pred_probs)
#' summary(accuracy(classified_LT5SR_1986$model))
#' }
classify_image <- function(x, train_data, pred_classes_filename=NULL, 
                           pred_probs_filename=NULL, classProbs=TRUE, 
                           use_training_flag=TRUE, tune_length=8, 
                           train_control=NULL, tune_grid=NULL) {

    cl <- options('rasterClusterObject')[[1]]
    if (is.null(cl)) {
        inparallel <- FALSE
    } else {
        if (!require(doSNOW)) {
            warning('Cluster object found, but "doSNOW" package is required to run training in parallel. Running sequentially.')
        }
        registerDoSNOW(cl)
        inparallel <- TRUE
    }

    message('Training classifier...')
    if (is.null(train_control)) {
        train_control <- trainControl(method="repeatedcv", repeats=5, 
                                      classProbs=classProbs)
    }
    if (use_training_flag) {
        if (!('Training' %in%names(train_data))) {
            stop('when use_training_flag is TRUE, train_data must have a "Training" column')
        }
    }
    # Build the formula, excluding the training flag column (if it exists) from 
    # the model formula
    formula_vars <- names(train_data)
    formula_vars <- formula_vars[!(formula_vars %in% c('y', 'Training', 
                                                       'Poly_FID'))]
    model_formula <- formula(paste('y ~', paste(formula_vars, collapse=' + ')))

    model <- train(model_formula, data=train_data, method="svmRadial",
                   preProc=c('center', 'scale'), subset=train_data$Training,
                   tuneLength=tune_length, trControl=train_control, 
                   tuneGrid=tune_grid)

    message('Predicting classes...')
    n_classes <- length(levels(model))
    if (inparallel) {
        pred_classes <- clusterR(x, predict, args=list(model))
    } else {
        pred_classes <- predict(x, model)
    }
    names(pred_classes) <- 'cover'

    if (classProbs) {
        message('Calculating class probabilities...')
        if (inparallel) {
            pred_probs <- clusterR(x, predict, args=list(model, type="prob", index=c(1:n_classes)))
        } else {
            pred_probs <- predict(x, model, type="prob", index=c(1:n_classes))
        }
        return(list(model=model, pred_classes=pred_classes, pred_probs=pred_probs))
    } else {
        return(list(model=model, pred_classes=pred_classes))
    }
}
