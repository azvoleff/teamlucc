#' Classifies an image using a trained random forest classifier
#'
#' @export
#' @import caret randomForest
#' @importFrom spatial.tools predict_rasterEngine
#' @param x a \code{Raster*} image with the predictor layer(s) for the 
#' classification
#' @param model a trained random forest model as output by 
#' \code{\link{rf_train}}
#' @param class_probs whether to also calculate and return the probabilities of 
#' membership for each class
#' @param use_training_flag indicates whether to exclude data flagged as 
#' testing data when training the classifier. For this to work the input 
#' train_data \code{data.frame} must have a column named 'training_flag' that 
#' indicates, for each pixel, whether that pixel is a training pixel (coded as 
#' TRUE) or testing pixel (coded as FALSE).
#' @param notify notifier to use (defaults to \code{print} function). See the 
#' \code{notifyR} package for one way of sending notifications from R. The 
#' \code{notify} function should accept a string as the only argument.
#' @return a list with 2 elements: the predicted classes as a 
#' \code{RasterLayer} and the class probabilities as a \code{RasterBrick}
#' @examples
#' #TODO: Add example.
rf_classify <- function(x, model, class_probs=TRUE, notify=print) {
    if (class_probs) {
        pred_probs <- predict_rasterEngine(object=model, newdata=x, 
                                           type='prob')
        names(pred_probs) <- levels(model)
        # TODO: Need to calculate the index of the highest probability 
        # prediction for each pixel and assign that value to the pred_classes 
        # layer
    } else {
        pred_classes <- predict_rasterEngine(object=model, newdata=x)
        pred_probs <- NULL
    }
    names(pred_classes) <- 'prediction'

    return(list(pred_classes=pred_classes, pred_probs=pred_probs))
}
