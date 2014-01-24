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
#' @param output_path the path to use for the output @param n_cpus the number 
#' of CPUs to use for processes that can run in parallel
#' @param overwrite whether to overwrite existing files (otherwise an error 
#' will be raised)
#' @param notify notifier to use (defaults to \code{print} function). See the 
#' \code{notifyR} package for one way of sending notifications from R. The 
#' \code{notify} function should accept a string as the only argument.
#' @examples
#' #TODO: Add example
team_classify <- function(predictors, train_data, output_path, n_cpus=1, 
                          overwrite=FALSE, notify=print) {
    timer <- Track_time(notify)

    timer <- start_timer(timer, label='Running team_classify')

    if (n_cpus > 1) beginCluster(n_cpus)

    timer <- start_timer(timer, label='Running classify_image')
    classification <- classify_image(predictors, train_data, notify=notify)
    timer <- stop_timer(timer, label='Running classify_image')
    classification$model
    plot(classification$pred_classes)

    # Perform accuracy assessment using an independent dataset:
    timer <- start_timer(timer, label='Running accuracy assessment')
    acc <- accuracy(classification$model, 
                    pop=classification$pred_classes)
    summary(acc)
    timer <- stop_timer(timer, label='Running accuracy assessment')

    if (n_cpus > 1) endCluster()

    timer <- stop_timer(timer, label='Running team_classify')
}
