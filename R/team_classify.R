#' Classify a preprocessed surface reflectance image
#'
#' First the image should be preprocesed using the \code{team_preprocess} 
#' function.
#'
#' @export
#' @param preds the path to an raster stack of predictor layers output by the 
#' \code{team_preprocess} function
team_classify <- function(preds) {
    train_polys_1986 <- readOGR('H:/Data/TEAM/VB/Vectors', 'VB_training_1986_037_LT5')
    #train_polys_2012 <- readOGR('H:/Data/TEAM/VB/Vectors', 'VB_training_2012_021_LE7')
    train_data_1986 <- extract_training_data(L5TSR_1986_preds, 
                                             train_polys_1986, "Poly_Type", 
                                             training=.7)
    classification_1986_L5TSR$model
    plot(classification_1986_L5TSR$pred_classes)

    # Perform accuracy assessment using an independent dataset:
    acc <- accuracy(classification_1986_L5TSR$model, 
                    pop=classification_1986_L5TSR$pred_classes)
    summary(acc)
}
