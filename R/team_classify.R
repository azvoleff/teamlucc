#' Classify a preprocessed surface reflectance image
#'
#' First the image should be preprocesed using the \code{team_preprocess} 
#' function.
#'
#' @export
#' @param predictors a \code{Raster*} of predictor layers output by the 
#' \code{team_preprocess} function or path to an image stack in a format 
#' readable by the \code{raster} package.
#' @param train_data a training dataset as output by 
#' \code{extract_training_data}
#' @param overwrite whether to overwrite existing files (otherwise an error 
#' will be raised)
#' @examples
#' #TODO: Add example
team_classify <- function(predictors, train_data, overwrite=FALSE) {
    classification <- classify_image(predictors, train_data)
    classification$model
    plot(classification$pred_classes)

    # Perform accuracy assessment using an independent dataset:
    acc <- accuracy(classification$model, 
                    pop=classification$pred_classes)
    summary(acc)
}
