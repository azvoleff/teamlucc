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
svm_classify <- function(x, train_data, out_file_base) {
    svm_tune <- tune.svm(y ~ ., data=train_data,
                         gamma=2^seq(-1, 1, .5), cost=2^seq(2, 4, .5),
                         probability=TRUE)
    svm_best <- svm_tune$best.model
    # Setup two rasters. out_classes wills store the predicted classes for the 
    # classified image. out_probs wills store a multiband raster with one band 
    # per class, with the predicted probabilities of class membership for each 
    # class.
    out_classes <- raster(x)
    out_classes <- writeStart(out_classes, paste(out_file_base, '_classes.envi', sep=''))
    out_probs <- brick(raster(x), values=FALSE, nl=svm_best$nclasses)
    out_probs <- writeStart(out_probs, paste(out_file_base, '_probs.envi', sep=''))
    # Process over blocks (rather than row by row) to save processing time.
    bs <- blockSize(x)
    for (block_num in 1:bs$n) {
        this_block <- data.frame(getValues(x, row=bs$row[block_num], nrows=bs$nrows[block_num]))
        names(this_block) <- names(x)
        # First write predicted classes
        pred <- predict(svm_best, newdata=this_block, probability=TRUE)
        writeValues(out_classes, as.numeric(pred), bs$row[block_num])
        # Now write predicted probabilities
        pred_probs <- attr(pred, 'prob')
        writeValues(out_probs, pred_probs, bs$row[block_num])
    }
    out_classes <- writeStop(out_classes)
    out_probs <- writeStop(out_probs)
    return(svm_best)
}
