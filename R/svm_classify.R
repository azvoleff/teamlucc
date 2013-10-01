#' Runs a Support Vector Machine (SVM) classification
#'
#' @export
#' @param x a \code{Raster*} image with the predictor layer(s) for the 
#' classification
#' @param train_data a data table with a column labeled 'y' with the observed 
#' classes, and one or more columns with the values of predictor(s) at each 
#' location.
#' @param out_file_base the base name to use when naming the output files
#' @return the best SVM chosen after tuning
#' @examples
#' L5TSR_1986 <- stack(system.file('extdata/L5TSR_1986.dat', package='teamr'))
#' data(L5TSR_1986_training)
#' train_data_1986 <- extract_training_data(L5TSR_1986, L5TSR_1986_training)
#' \dontrun{
#' # This code is not run because it writes to a local folder
#' svm_classify(L5TSR_1986, train_data_1986, 'L5TSR_1986_classified')
#' }
svm_classify <- function(x, train_data, out_file_base,                          
                         pred_classes_filename=NULL, pred_prob_filename=NULL, 
                         train_grid=NULL) {
    print('Training SVM...')
    if (is.null(train_grid)) {
        sig_dist <- as.vector(sigest(y ~ ., data=train_data, frac=1))
        svm_train_grid <- data.frame(.sigma=sig_dist[1], .C=2^(-2:12))
    }
    svm_train_control <- trainControl(method="repeatedcv",
                                      repeats=5,
                                      classProbs=TRUE)
    svm_train <-  train(y ~ ., data=train_data, method="svmRadial",
                       preProc=c('center', 'scale'),
                       tuneGrid=svm_train_grid, trControl=svm_train_control)

    print('Predicting classes...')
    calc_pred_classes <- function(x, svm_train, ...) {
        preds <- predict(x, svm_train)
        preds <- array(getValues(preds), dim=c(dim(x)[1], dim(x)[2], 1))
        return(preds)
    }
    pred_classes <- focal_hpc(x=x, fun=calc_pred_classes, args=list(svm_train), 
                     filename=pred_classes_filename, chunk_format="raster")

    print('Predicting class probabilities...')
    calc_pred_probs <- function(x, svm_train, ...) {
        preds <- predict(x, svm_train, type="prob", index=c(1:dim(x)[3]))
        preds <- array(getValues(preds), dim=c(dim(x)[1], dim(x)[2], dim(x)[3]))
        return(preds)
    }
    pred_probs <- focal_hpc(x=x, fun=calc_pred_probs, args=list(svm_train), 
                            filename=pred_prob_filename, chunk_format="raster")

    return(list(svm_train=svm_train, pred_classes=pred_classes, pred_probs=pred_probs))
}
