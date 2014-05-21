#' Runs an image classification using an SVM classifier
#'
#' Currently only supports classification using a support vector machine (SVM).
#'
#' Processing can be done in parallel using all using the cluster facilities in 
#' the \code{spatial.tools} package. To enable clustering, call 
#' \code{beginCluster} before running \code{classify_image}.  To stop the 
#' cluster when finished, call \code{endCluster}.
#'
#' @export
#' @import caret e1071 kernlab randomForest
#' @param x a \code{Raster*} image with the predictor layer(s) for the 
#' classification
#' @param train_data a \code{link{pixel_data}} object
#' @param class_probs whether to also calculate and return the probabilities of 
#' membership for each class
#' @param use_training_flag indicates whether to exclude data flagged as 
#' testing data when training the classifier. For this to work the input 
#' train_data \code{data.frame} must have a column named 'training_flag' that 
#' indicates, for each pixel, whether that pixel is a training pixel (coded as 
#' TRUE) or testing pixel (coded as FALSE).
#' @param tune_length the number of levels of each parameter that should be 
#' tried by \code{train} when training the classifier
#' @param train_control default is NULL (reasonable value will be set 
#' automatically).  For details see \code{\link{trainControl}}.
#' @param tune_grid the training grid to be used for training the classifier.  
#' Must be a \code{data.frame} with two columns: ".sigma" and ".C".
#' @param do_split whether to use normal mixture modeling to split the 
#' input classes into subsets to aid in classifying spectrally diverse classes
#' @param use_rfe whether to use Recursive Feature Extraction (RFE) as 
#' implemented in the \code{caret} package to select a subset of the input 
#' features to be used in the classification
#' @param notify notifier to use (defaults to \code{print} function). See the 
#' \code{notifyR} package for one way of sending notifications from R. The 
#' \code{notify} function should accept a string as the only argument.
#' @return a list with 3 elements: the trained classifier, the predicted classes 
#' \code{RasterLayer} and the class probabilities \code{RasterBrick}
#' @examples
#' \dontrun{
#' # Don't run long example
#' set.seed(0)
#' train_data_1986 <- get_pixels(L5TSR_1986, polys=L5TSR_1986_2001_training,
#'                               class_col="class_1986", training=.7)
#' L5TSR_1986_classified <- classify_image(L5TSR_1986, train_data_1986)
#'
#' L5TSR_1986_classified$model
#' plot(L5TSR_1986_classified$pred_classes)
#' plot(L5TSR_1986_classified$pred_probs)
#' accuracy(L5TSR_1986_classified$model)
#' }
classify_image <- function(x, train_data, class_probs=TRUE, 
                           use_training_flag=TRUE, tune_length=8, 
                           train_control=NULL, tune_grid=NULL, 
                           do_split=FALSE, use_rfe=FALSE, notify=print) {
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
        train_control <- trainControl(method="repeatedcv", repeats=5, 
                                      classProbs=class_probs)
    }
    # Build the formula, excluding the training flag column (if it exists) from 
    # the model formula
    predictor_names <- names(train_data@x)
    model_formula <- formula(paste('y ~',
                                   paste(predictor_names, collapse=' + ')))

    if (do_split) {
        notify('Performing pre-classification clustering...')
        training_split <- split_classes(train_data)
        y <- training_split$y

    } else {
        training_split <- NULL
        y <- train_data@y
    }

    if (use_rfe) {
        stop('rfe not yet supported')
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

    train_data <- cbind(y=y,
                        train_data@x,
                        training_flag=train_data@training_flag,
                        poly_ID=train_data@poly_ID)

    notify('Training classifier...')
    model <- train(model_formula, data=train_data, method="svmRadial",
                   preProc=c('range'), subset=train_data$training_flag,
                   tuneLength=tune_length, trControl=train_control, 
                   tuneGrid=tune_grid)

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

    if (do_split) {
        notify('Recoding split classes...')
        is_becomes <- cbind(training_split$reclass_mat$split_id,
                            training_split$reclass_mat$id)
        pred_classes_recode <- reclassify(pred_classes, is_becomes)
        pred_classes_recode <- ratify(pred_classes_recode)
        rat <- levels(pred_classes_recode)[[1]]
        rat$value <- unique(training_split$reclass_mat$name)
        levels(pred_classes_recode) <- rat
        names(pred_classes_recode) <- names(pred_classes)
        if (class_probs) {
            pred_probs_recode <- stackApply(pred_probs, is_becomes[, 2], sum)
            names(pred_probs_recode) <- unique(training_split$reclass_mat$name)
        }
    } else {
        pred_classes_recode <- NULL
        pred_probs_recode <- NULL
    }

    return(list(model=model, pred_classes=pred_classes, pred_probs=pred_probs, 
                split_classes=training_split, 
                pred_probs_recode=pred_probs_recode, 
                pred_classes_recode=pred_classes_recode))
}
