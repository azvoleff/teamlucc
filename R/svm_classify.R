#' Runs a Support Vector Machine (SVM) classification
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
#' @param train_grid the training grid to be used for training the SVM. Must be 
#' a \code{data.frame} with two columns: ".sigma" and ".C".
#' @param classProbs whether to also calculate and return the probabilities of 
#' membership for each class
#' @param use_training_flag indicates whether to exclude data flagged as 
#' testing data when training the classifier. For this to work the input 
#' train_data \code{data.frame} must have a column named 'Training' that 
#' indicates, for each pixel, whether that pixel is a training pixel (coded as 
#' TRUE) or testing pixel (coded as FALSE).
#' @return a list with 3 elements: the trained SVM model, the predicted classes 
#' \code{RasterLayer} and the class probabilities \code{RasterBrick}
#' @details processing can be done in parallel using all available CPUs through 
#' the use of the cluster facilities in the \code{spatial.tools} package. To 
#' enable clustering, call \code{sfQuickInit} before running 
#' \code{svm_classify}. To stop the cluster when finished, call 
#' \code{sfQuickStop}.
#'
#' @examples
#' \dontrun{
#' # Don't run long example
#' # Load example datasets
#' L5TSR_1986 <- stack(system.file('extdata/L5TSR_1986.dat', package='teamr'))
#' train_data_1986 <- extract_training_data(L5TSR_1986, 
#'                                          polys=L5TSR_1986_2001_training,
#'                                          classcol="t1_class", testfrac=.3)
#' # Supply a small training grid for svm_classify to save processing time for 
#' # the purposes of this example - in normal use, train_grid can be left 
#' # unspecified.
#' svmout <- svm_classify(L5TSR_1986, train_data_1986, 
#'                        train_grid=data.frame(.sigma=.0495, .C=0.5),
#'                        classProbs=TRUE)
#'
#' # Examine output from svm_classify
#' svmout$svm_train
#' plot(svmout$pred_classes)
#' plot(svmout$pred_probs)
#' }
svm_classify <- function(x, train_data, pred_classes_filename=NULL, 
                         pred_probs_filename=NULL, train_grid=NULL,
                         classProbs=FALSE, use_training_flag=TRUE) {
    message('Training SVM...')
    if (is.null(train_grid)) {
        sig_dist <- as.vector(sigest(y ~ ., data=train_data, frac=1))
        train_grid <- data.frame(.sigma=sig_dist[1], .C=2^(-6:12))
    }
    if (use_training_flag) {
        if (!('Training' %in%names(train_data))) {
            stop('when use_training_flag is TRUE, train_data must have a "Training" column')
        }
        train_data <- train_data[train_data$Training, ]
    }
    svm_train_control <- trainControl(method="repeatedcv",
                                      repeats=5,
                                      classProbs=classProbs)
    # Build the formula, excluding the training flag column (if it exists) from 
    # the model formula
    formula_vars <- names(train_data)
    formula_vars <- formula_vars[!(formula_vars %in% c('y', 'Training'))]
    svm_formula <- formula(paste('y ~', paste(formula_vars, collapse=' + ')))
    svm_train <-  train(svm_formula, data=train_data, method="svmRadial",
                        preProc=c('center', 'scale'), subset=train_data_1986$Training,
                        tuneGrid=train_grid, trControl=svm_train_control)

    message('Predicting classes...')
    calc_preds <- function(in_rast, svm_train, n_classes, classProbs, ...) {
        preds <- predict(in_rast, svm_train)
        if (classProbs) {
            pred_probs <- predict(in_rast, svm_train, type="prob", 
                                  index=c(1:n_classes))
            preds <- stack(preds, pred_probs)
        }
        preds <- array(getValues(preds), dim=c(dim(in_rast)[2], 
                                               dim(in_rast)[1], 
                                               nlayers(preds)))
        return(preds)
    }
    n_classes <- length(caret:::getClassLevels(svm_train))
    preds <- rasterEngine(in_rast=x, fun=calc_preds, 
                               args=list(svm_train=svm_train, 
                                         n_classes=n_classes,
                                         classProbs=classProbs), 
                               filename=pred_probs_filename, 
                               chunk_format="raster")
    pred_classes <- raster(preds, layer=1)
    names(pred_classes) <- 'cover'
    if (classProbs) {
        pred_probs <- dropLayer(preds, 1)
        names(pred_probs) <- caret:::getClassLevels(svm_train)
        return(list(svm_train=svm_train, pred_classes=pred_classes, pred_probs=pred_probs))
    } else {
        return(list(svm_train=svm_train, pred_classes=pred_classes))
    }

}
