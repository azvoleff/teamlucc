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
    # Assign standardized layer names to input image so that different images 
    # can be used with the same model
    names(x) <- paste0('pred', seq(1:nlayers(x)))

    notify('Predicting classes...')
    n_classes <- length(levels(model))
    if (inparallel) {
        pred_classes <- clusterR(x, predict, args=list(model))
    } else {
        pred_classes <- predict_rasterEngine(object=model, newdata=x)
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

    return(list(pred_classes=pred_classes, pred_probs=pred_probs))
}
