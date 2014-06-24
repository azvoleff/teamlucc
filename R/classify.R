#' Classifies an image using a trained random forest classifier
#'
#' @export
#' @import caret randomForest
#' @importFrom spatial.tools rasterEngine predict_rasterEngine
#' @param x a \code{Raster*} image with the predictor layer(s) for the 
#' classification
#' @param model a trained random forest model as output by 
#' \code{\link{rf_train}}
#' @param classes_file filename for predicted classes (or missing)
#' @param prob_file filename for predicted probabilities (or missing) CURRENTLY 
#' IGNORED
#' @param use_training_flag indicates whether to exclude data flagged as 
#' testing data when training the classifier. For this to work the input 
#' train_data \code{data.frame} must have a column named 'training_flag' that 
#' indicates, for each pixel, whether that pixel is a training pixel (coded as 
#' TRUE) or testing pixel (coded as FALSE).
#' @param overwrite whether to overwrite \code{out_name} if it already exists
#' @param notify notifier to use (defaults to \code{print} function). See the 
#' \code{notifyR} package for one way of sending notifications from R. The 
#' \code{notify} function should accept a string as the only argument.
#' @return a list with 2 elements: the predicted classes as a 
#' \code{RasterLayer} and the class probabilities as a \code{RasterBrick}
#' @examples
#' #TODO: Fix example
#' model <- rf_train(x, tr_data)
#' preds <- rf_classify(x, model)
#' require(rasterVis)
#' levelplot(preds$classes)
classify <- function(x, model, classes_file, prob_file, overwrite=FALSE, 
                     notify=print) {
    #if (missing(prob_file)) prob_file <- rasterTmpFile()
    #TODO: find out how to use prob_file with spatial.tools
    probs <- predict_rasterEngine(object=model, newdata=x, type='prob')
    names(probs) <- levels(model)

    # Calculate the highest probability class from the class probabilities
    if (missing(classes_file)) classes_file <- rasterTmpFile()
    classes <- calc(probs, fun=function(vals) {
        # Subtract 1 below as software like ENVI starts class codes at zero
        out <- as.numeric(which(vals == max(vals))) - 1
        #TODO: Need to handle case of ties (length(out) > 1)
        if (length(out) != 1) out <- NA
        return(out)
    }, datatype='INT2S', filename=classes_file, overwrite=overwrite)
    names(classes) <- 'prediction'

    codes <- data.frame(code=seq(0, (nlevels(model) - 1)), class=levels(model))

    return(list(classes=classes, probs=probs, codes=codes))
}
