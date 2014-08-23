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
#' @importFrom spatial.tools rasterEngine
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
    # TODO: Check with Jonathan why below fix is needed
    if (!("RasterBrick" %in% class(x))) x <- brick(x)

    if (!missing(prob_file) && file_test('-f', prob_file) && !overwrite) {
        stop(paste('output file', prob_file, 'already exists and overwrite=FALSE'))
    }
    if (!missing(classes_file) && file_test('-f', classes_file) && !overwrite) {
        stop(paste('output file', classes_file, 'already exists and overwrite=FALSE'))
    }

    make_preds <- function(inrast, model, factors, ...) {
        # First, preserve the names:
        band_names <- dimnames(inrast)[3][[1]]

        # Flatten the array to a matrix (we lose the names here)
        inrast_mat <- inrast
        dim(inrast_mat) <- c(dim(inrast)[1]*dim(inrast)[2], dim(inrast)[3])
        inrast_df <- as.data.frame(inrast_mat)
        names(inrast_df) <- band_names

        # Make sure any factor variables are converted to factors and that the 
        # proper levels are assigned
        if (length(factors) > 0) {
            for (n in 1:length(factors)) {
                factor_var <- names(factors)[n]
                factor_col <- which(names(inrast_df) == factor_var)
                inrast_df[, factor_col] <- factor(inrast_df[, factor_col], 
                                                  levels=factors[[n]])
            }
        }

        good_obs <- complete.cases(inrast_df)
        preds <- matrix(NA, nrow=nrow(inrast_df), ncol=nlevels(model))
        if (sum(good_obs) > 0) {
            good_preds <- predict(model, inrast_df[good_obs, ], type="prob")
            preds[which(good_obs), ] <- as.matrix(good_preds)
        }

        preds_array <- array(preds, dim=c(dim(inrast)[1], dim(inrast)[2], 
                                          nlevels(model)))
        return(preds_array)
    }
    probs <- rasterEngine(inrast=x, fun=make_preds,
                          args=list(model=model, factors=factors),
                          filename=rasterTmpFile(), overwrite=overwrite, 
                          datatype="FLT4S", .packages=c("randomForest"),
                          setMinMax=TRUE)
    # spatial.tools can only output the raster package grid format - so output 
    # to a tempfile in that format then copy over to the requested final output 
    # format if a filename was supplied
    if (!missing(prob_file)) {
        probs <- writeRaster(probs, filename=prob_file, overwrite=overwrite, 
                             datatype='FLT4S')
    }
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
