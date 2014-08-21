#' Classify an image using a trained classifier
#'
#' This function will produce two outputs - a prediction image and a 
#' probability image. The prediction image contains the predicted classes, the 
#' and the probability image contains the per-pixel predicted probabilities of 
#' occurrence of each class.
#'
#' This function will run in parallel if a parallel backend is registered with 
#' \code{\link{foreach}} - TEMPORARILY DISABLED.
#'
#' @export
#' @import caret
#' @importFrom spatial.tools rasterEngine predict_rasterEngine
#' @param x a \code{Raster*} image with the predictor layer(s) for the 
#' classification
#' @param model a trained classifier as output by 
#' \code{\link{train_classifier}}
#' @param classes_file filename for predicted classes (or missing)
#' @param prob_file filename for predicted probabilities (or missing)
#' @param factors a list of character vector giving the names of predictors 
#' (layer names from the images used to build \code{train_data}) that should be 
#' treated as factors, and specifying the levels of each factor. For example, 
#' \code{factors=list(year=c(1990, 1995, 2000, 2005, 2010))}.
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
classify <- function(x, model, classes_file, prob_file, factors=list(), 
                     overwrite=FALSE) {
    if (!missing(prob_file) && file_test('-f', prob_file) && !overwrite) {
        stop(paste('output file', prob_file, 'already exists and overwrite=FALSE'))
    }
    if (!missing(classes_file) && file_test('-f', classes_file) && !overwrite) {
        stop(paste('output file', classes_file, 'already exists and overwrite=FALSE'))
    }

    # # Below is currently disabled due to issues handling factors in
    # # predict_rasterEngine
    # probs <- predict_rasterEngine(object=model, newdata=x, type='prob')
    #
    # # spatial.tools can only output the raster package grid format - so output 
    # # to a tempfile in that format then copy over to the requested final output 
    # # format if a filename was supplied
    # if (!missing(prob_file)) {
    #     probs <- writeRaster(probs, filename=prob_file, overwrite=overwrite, 
    #                          datatype='FLT4S')
    # }

    if (missing(prob_file)) {
        prob_file <- rasterTmpFile()
    }
    probs <- predict(x, model, type="prob", index=c(1:nlevels(model)), 
                     factors=factors, filename=prob_file, overwrite=overwrite, 
                     datatype="FLT4S")

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
