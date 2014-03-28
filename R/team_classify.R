#' Classify a preprocessed surface reflectance image
#'
#' First the image should be preprocesed using the \code{team_preprocess} 
#' function.
#'
#' @export
#' @importFrom rgdal readOGR
#' @importFrom sp spTransform
#' @importFrom tools file_path_sans_ext
#' @param predictor_file a \code{Raster*} of predictor layers output by the 
#' \code{team_preprocess} function or path to an image stack in a format 
#' readable by the \code{raster} package.
#' @param train_shp a training dataset as output by 
#' \code{extract_observed}
#' @param output_path the path to use for the output
#' @param class_col the name of the column containing the response variable 
#' (for example the land cover type of each pixel)
#' @param training indicator of which polygons to use in training. Can be: 1) a 
#' string giving the name of a column indicating whether each polygon is to be 
#' used in training (column equal to TRUE) or in testing (column equal to 
#' FALSE), or 2) a logical vector of length equal to length(polys), or 3) a 
#' number between 0 and 1 indicating the fraction of the polygons to be 
#' randomly selected for use in training.
#' @param n_cpus the number of CPUs to use for processes that can run in 
#' parallel
#' @param overwrite whether to overwrite existing files (otherwise an error 
#' will be raised)
#' @param notify notifier to use (defaults to \code{print} function). See the 
#' \code{notifyR} package for one way of sending notifications from R. The 
#' \code{notify} function should accept a string as the only argument.
#' @examples
#' #TODO: Add example
team_classify <- function(predictor_file, train_shp, output_path, 
                          class_col="Poly_Type", training=.6, n_cpus=1, 
                          overwrite=FALSE, notify=print) {
    if (!file_test("-f", train_shp)) {
        stop(paste(train_shp, "does not exist"))
    }
    if (!file_test("-f", predictor_file)) {
        stop(paste(predictor_file, "does not exist"))
    }
    if (!file_test("-d", output_path)) {
        stop(paste(output_path, "does not exist"))
    }

    if (n_cpus > 1) beginCluster(n_cpus)

    timer <- Track_time(notify)
    timer <- start_timer(timer, label='Running team_classify')

    predictors <- brick(predictor_file)
    pred_rast_basename <- basename(file_path_sans_ext(predictor_file))

    train_polys <- readOGR(dirname(train_shp), basename(file_path_sans_ext(train_shp)))
    train_polys <- spTransform(train_polys, crs(predictors))
    train_data <- extract_observed(predictors, train_polys, 
                                        class_col=class_col, training=training)

    timer <- start_timer(timer, label='Running classify_image')
    classification <- classify_image(predictors, train_data, notify=notify)
    model <- classification$model
    save(model, file=file.path(output_path, paste(pred_rast_basename, 
                                                  'predmodel.RData', sep='_')))
    writeRaster(classification$pred_classes,
                filename=file.path(output_path, paste(pred_rast_basename, 
                                                      'predclasses.envi', 
                                                      sep='_')),
                datatype='INT2S', overwrite=overwrite)
    writeRaster(scale_raster(classification$pred_probs),
                filename=file.path(output_path, paste(pred_rast_basename, 
                                                      'predprobs.envi', 
                                                      sep='_')),
                datatype='INT2S', overwrite=overwrite)
    timer <- stop_timer(timer, label='Running classify_image')

    # cls <- levels(train_data$y) 
    # cls <- data.frame(code=seq(1:length(cls)), name=cls)
    # color_image(classification$predclasses, cls,
    #             file.path(output_path, paste(pred_rast_basename, 
    #             'predclasses_colored.envi', sep='_')))

    # Perform accuracy assessment using an independent dataset:
    timer <- start_timer(timer, label='Running accuracy assessment')
    acc <- accuracy(classification$model, 
                    pop=classification$pred_classes)
    capture.output(summary(acc),
                   file=file.path(output_path, paste(pred_rast_basename, 'predacc.txt', sep='_')))
    timer <- stop_timer(timer, label='Running accuracy assessment')

    if (n_cpus > 1) endCluster()

    timer <- stop_timer(timer, label='Running team_classify')
}
