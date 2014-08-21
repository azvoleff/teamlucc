#' Classify an image using a trained classifier
#'
#' This function will produce two outputs - a prediction image and a 
#' probability image. The prediction image contains the predicted classes, the 
#' and the probability image contains the per-pixel predicted probabilities of 
#' occurrence of each class.
#'
#' This function will run in parallel if a parallel backend is registered with 
#' \code{\link{foreach}}.
#'
#' @export
#' @import caret
#' @importFrom spatial.tools rasterEngine predict_rasterEngine
#' @param x a \code{Raster*} image with the predictor layer(s) for the 
#' classification
#' @param model a trained classifier as output by 
#' \code{\link{train_classifier}}
#' @param classes_file filename for predicted classes (or missing)
#' @param prob_file filename for predicted probabilities (or missing) CURRENTLY 
#' IGNORED
#' @param overwrite whether to overwrite \code{out_name} if it already exists
#' @return a list with 2 elements: the predicted classes as a 
#' \code{RasterLayer} and the class probabilities as a \code{RasterBrick}
#' @examples
#' \dontrun{
#' train_data <- get_pixels(L5TSR_1986, L5TSR_1986_2001_training, "class_1986", 
#'                          training=.6)
#' model <- train_classifier(train_data)
#' preds <- classify(L5TSR_1986, model)
#' plot(preds$classes)
#' plot(preds$probs)
#' }
classify <- function(x, model, classes_file, prob_file, overwrite=FALSE) {
    if (!missing(prob_file) && file_test('-f', prob_file) && !overwrite) {
        stop(paste('output file', prob_file, 'already exists and overwrite=FALSE'))
    }
    if (!missing(classes_file) && file_test('-f', classes_file) && !overwrite) {
        stop(paste('output file', classes_file, 'already exists and overwrite=FALSE'))
    }

    # probs <- predict_rasterEngine(object=model, newdata=x, type='prob')
    #
    # spatial.tools can only output the raster package grid format - so output 
    # to a tempfile in that format then copy over to the requested final output 
    # format if a filename was supplied
    # if (!missing(prob_file)) {
    #     probs <- writeRaster(probs, filename=prob_file, overwrite=overwrite, 
    #                          datatype='FLT4S')
    # }

    if (missing(prob_file)) {
        prob_file <- rasterTmpFile()
    }
    probs <- predict(x, model,  type="prob", progress='text', 
                     index=c(1:nlevels(model)), filename=prob_file, 
                     overwrite=overwrite, datatype="FLT4S")

    # probs <- predict(x, model,  type="prob", progress='text',
    #                  factors=list(aspect=c(1, 2, 3, 4),
    #                               year=c(1990, 1995, 2000, 2005, 2010)),
    #                  index=c(1:nlevels(model)), filename=prob_file, 
    #                  overwrite=overwrite, datatype="FLT4S")

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
