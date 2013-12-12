#' Runs an image classification
#'
#' Currently only supports classification using a support vector machine (SVM).
#'
#' @export
#' @import kernlab e1071 caret spatial.tools
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
#' @param train_grid the training grid to be used for training the classifier. Must be 
#' a \code{data.frame} with two columns: ".sigma" and ".C".
#' @param classProbs whether to also calculate and return the probabilities of 
#' membership for each class
#' @param use_training_flag indicates whether to exclude data flagged as 
#' testing data when training the classifier. For this to work the input 
#' train_data \code{data.frame} must have a column named 'Training' that 
#' indicates, for each pixel, whether that pixel is a training pixel (coded as 
#' TRUE) or testing pixel (coded as FALSE).
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
#' train_data_1986 <- extract_training_data(L5TSR_1986, 
#'                                          polys=L5TSR_1986_2001_training,
#'                                          classcol="t1_class", training=.7)
#' # Supply a small training grid for classify_image to save processing time for 
#' # the purposes of this example - in normal use, train_grid can be left 
#' # unspecified.
#' classified_LT5SR_1986 <- classify_image(L5TSR_1986, train_data_1986, 
#'     train_grid=data.frame(.sigma=.0495, .C=0.5), classProbs=TRUE)
#'
#' classified_LT5SR_1986$model
#' plot(classified_LT5SR_1986$pred_classes)
#' plot(classified_LT5SR_1986$pred_probs)
#' }
classify_image <- function(x, train_data, pred_classes_filename=NULL, 
                           pred_probs_filename=NULL, train_grid=NULL,
                           classProbs=FALSE, use_training_flag=TRUE) {
    message('Training classifier...')
    if (is.null(train_grid)) {
        sig_dist <- as.vector(sigest(y ~ ., data=train_data, frac=1))
        train_grid <- data.frame(.sigma=sig_dist[1], .C=2^(-6:12))
    }
    if (use_training_flag) {
        if (!('Training' %in%names(train_data))) {
            stop('when use_training_flag is TRUE, train_data must have a "Training" column')
        }
    }
    svm_train_control <- trainControl(method="repeatedcv",
                                      repeats=5,
                                      classProbs=classProbs)
    # Build the formula, excluding the training flag column (if it exists) from 
    # the model formula
    formula_vars <- names(train_data)
    formula_vars <- formula_vars[!(formula_vars %in% c('y', 'Training', 'Poly_FID'))]
    model_formula <- formula(paste('y ~', paste(formula_vars, collapse=' + ')))
    model <- train(model_formula, data=train_data, method="svmRadial",
                   preProc=c('center', 'scale'), subset=train_data$Training,
                   tuneGrid=train_grid, trControl=svm_train_control)

    message('Predicting classes...')
    calc_preds <- function(in_rast, model, n_classes, classProbs, ...) {
        preds <- predict(in_rast, model)
        if (classProbs) {
            pred_probs <- predict(in_rast, model, type="prob", 
                                  index=c(1:n_classes))
            preds <- stack(preds, pred_probs)
        }
        preds <- array(getValues(preds), dim=c(dim(in_rast)[2], 
                                               dim(in_rast)[1], 
                                               nlayers(preds)))
        return(preds)
    }
    n_classes <- length(caret:::getClassLevels(model))
    preds <- rasterEngine(in_rast=x, fun=calc_preds, 
                               args=list(model=model, 
                                         n_classes=n_classes,
                                         classProbs=classProbs), 
                               filename=pred_probs_filename, 
                               chunk_format="raster")
    pred_classes <- raster(preds, layer=1)
    names(pred_classes) <- 'cover'
    if (classProbs) {
        pred_probs <- dropLayer(preds, 1)
        names(pred_probs) <- caret:::getClassLevels(model)
        return(list(model=model, pred_classes=pred_classes, pred_probs=pred_probs))
    } else {
        return(list(model=model, pred_classes=pred_classes))
    }
}
