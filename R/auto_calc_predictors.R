#' Calculate predictor layers for a classification
#'
#' This function automates the calculation of a layer stack of predictor layers 
#' to use in a land use and/or land cover classification. See Details for the 
#' output layers.
#' 
#' The layers in the output layer stack are listed below. Note that all the 
#' layers are rescaled so that they range between -32,767 and 32,767 (allowing 
#' them to be stored as 16 bit unsigned integers).
#'
#' \bold{Predictor layer stack:}
#' \tabular{ll}{
#'     Layer 1: \tab Band 1 reflectance \cr
#'     Layer 2: \tab Band 2 reflectance \cr
#'     Layer 3: \tab Band 3 reflectance \cr
#'     Layer 4: \tab Band 4 reflectance \cr
#'     Layer 5: \tab Band 5 reflectance \cr
#'     Layer 6: \tab Band 7 reflectance \cr
#'     Layer 7: \tab MSAVI2 \cr
#'     Layer 8: \tab GLCM mean (from MSAVI2) \cr
#'     Layer 9: \tab GLCM variance (from MSAVI2) \cr
#'     Layer 10: \tab GLCM dissimilarity (from MSAVI2) \cr
#'     Layer 11: \tab Elevation \cr
#'     Layer 12: \tab Slope (radians X 10000) \cr
#'     Layer 13: \tab Aspect (see below) \cr
#' }
#'
#' The aspect is recoded as:
#'
#' \bold{Aspect coding:}
#' \tabular{ll}{
#'     1: \tab north facing (0-45 degrees, 315-360 degrees) \cr
#'     2: \tab east facing (45-135 degrees) \cr
#'     3: \tab south facing (135-225 degrees) \cr
#'     4: \tab west facing (225-315 degrees) \cr
#' }
#' @export
#' @importFrom glcm glcm
#' @importFrom stringr str_extract
#' @param x path to a preprocessed image as output by 
#' \code{auto_preprocess_landsat} or \code{auto_cloud_fill}.
#' @param dem DEM \code{RasterLayer} as output by \code{auto_setup_dem}
#' @param slopeaspect \code{RasterStack} as output by \code{auto_setup_dem}
#' @param output_path the path to use for the output (optional - if NULL then 
#' output images will be saved alongside the input images in the same folder).
#' @param ext file extension to use when saving output rasters (determines 
#' output file format).
#' @param overwrite whether to overwrite existing files (otherwise an error 
#' will be raised)
#' @param ...  additional arguments passed to \code{\link{glcm}}, such as
#' \code{n_grey}, \code{window}, or \code{shift}
#' @param notify notifier to use (defaults to \code{print} function). See the 
#' \code{notifyR} package for one way of sending notifications from R. The 
#' \code{notify} function should accept a string as the only argument.
auto_calc_predictors <- function(x, dem, slopeaspect, output_path=NULL, 
                                 ext='tif', overwrite=FALSE, notify=print,
                                 ...) {
    if (!file_test("-f", x)) {
        stop(paste("input image", x, "does not exist"))
    }
    if (!is.null(output_path) && !file_test("-d", output_path)) {
        stop(paste(output_path, "does not exist"))
    }

    ext <- gsub('^[.]', '', ext)

    timer <- Track_time(notify)
    timer <- start_timer(timer, label='Predictor calculation')

    # Setup a regex to identify preprocessed images
    preproc_regex <- '^[a-zA-Z]{2,3}_[0-9]{3}-[0-9]{3}_[0-9]{4}-[0-9]{3}_L[457][ET]SR(_tc)?'

    # Update the basename to refer the chosen file
    image_basename <- basename(file_path_sans_ext(x))

    image_stack <- brick(x)

    if (is.null(output_path)) {
        output_path <- dirname(x)
    }

    mask_stack_file <- paste0(file_path_sans_ext(x), '_masks.', ext)
    if (!file_test('-f', mask_stack_file)) {
        mask_stack_file <- gsub(paste0('(_tc)?.', ext, '$'), paste0('_masks.', ext), x)
        if (file_test('-f', mask_stack_file)) {
            warning('using masks file with old format (pre v0.5) teamlucc naming')
        } else {
            stop('could not find masks file')
        }
    }
    mask_stack <- brick(mask_stack_file)
    image_mask <- calc(mask_stack[[2]], function(maskvals) {
        # Mask clouds, cloud shadow, and fill
        (maskvals == 2) | (maskvals == 4) | (maskvals == 255)
    })

    ######################################################################
    # Calculate additional predictor layers (MSAVI and textures)
    timer <- start_timer(timer, label='Calculating MSAVI2')
    MSAVI2_filename <- file.path(output_path,
                                 paste0(image_basename, '_MSAVI2.', ext))
    MSAVI2_layer <- MSAVI2(red=raster(image_stack, layer=3),
                           nir=raster(image_stack, layer=4))
    # Truncate MSAVI2 to range between 0 and 1, and scale by 10,000 so it 
    # can be saved as a INT2S
    MSAVI2_layer <- calc(MSAVI2_layer, fun=function(vals) {
            vals[vals > 1] <- 1
            vals[vals < 0] <- 0
            vals <- round(vals * 10000)
        }, filename=MSAVI2_filename, overwrite=overwrite, datatype="INT2S")
    timer <- stop_timer(timer, label='Calculating MSAVI2')

    timer <- start_timer(timer, label='Calculating GLCM textures')
    MSAVI2_glcm_filename <- file.path(output_path,
                                      paste0(image_basename, 
                                            '_MSAVI2_glcm.', ext))
    glcm_statistics <- c('mean', 'variance', 'homogeneity', 'contrast', 
                         'dissimilarity', 'entropy', 'second_moment', 
                         'correlation')
    MSAVI2_layer[image_mask] <- NA
    # Need to know window and shift to calculate edge for apply_windowed. So if 
    # they are not in the dotted args, assume the defaults (since glcm will use 
    # the defaults if these parameters are not supplied).
    dots <- list(...)
    if (!("window" %in% names(dots))) {
        dots$window <- c(3, 3)
    }
    if (!("shift" %in% names(dots))) {
        dots$shift <- c(1, 1)
    }
    edge <- calc_glcm_edge(dots$shift, dots$window)
    # Note the min_x and max_x are given for MSAVI2 that has been scaled by 
    # 10,000
    apply_windowed_args <- list(x=MSAVI2_layer, fun=glcm, edge=edge, min_x=0, 
                             max_x=10000, filename=MSAVI2_glcm_filename, 
                             overwrite=overwrite, statistics=glcm_statistics, 
                             na_opt='center')
    apply_windowed_args <- c(apply_windowed_args, dots)
    MSAVI2_glcm <- do.call(apply_windowed, apply_windowed_args)
    names(MSAVI2_glcm) <- paste('glcm', glcm_statistics, sep='_')
    timer <- stop_timer(timer, label='Calculating GLCM textures')

    if (!missing(slopeaspect)) {
        timer <- start_timer(timer, label='Processing slopeaspect')
        names(slopeaspect) <- c('slope', 'aspect')
        # Classify aspect into north facing, east facing, etc., recalling 
        # that the aspect is stored in radians scaled by 1000.
        #     1: north facing (0-45, 315-360)
        #     2: east facing (45-135)
        #     3: south facing (135-225)
        #     4: west facing (225-315)
        aspect_cut <- raster::cut(slopeaspect$aspect/1000,
                                  c(-1, 45, 135, 225, 315, 361)*(pi/180))
        # Code both 0-45 and 315-360 aspect as North facing (1)
        aspect_cut[aspect_cut == 5] <- 1
        names(aspect_cut) <- 'aspect'
        timer <- stop_timer(timer, label='Processing slopeaspect')
    }

    ######################################################################
    # Layer stack predictor layers:
    timer <- start_timer(timer, label='Writing predictors')
    predictors <- stack(raster(image_stack, layer=1),
                        raster(image_stack, layer=2),
                        raster(image_stack, layer=3),
                        raster(image_stack, layer=4),
                        raster(image_stack, layer=5),
                        raster(image_stack, layer=6),
                        MSAVI2_layer,
                        scale_raster(MSAVI2_glcm$glcm_mean),
                        scale_raster(MSAVI2_glcm$glcm_variance),
                        scale_raster(MSAVI2_glcm$glcm_dissimilarity))
    predictor_names <- c('b1', 'b2', 'b3', 'b4', 'b5', 'b7', 'msavi', 
                         'msavi_glcm_mean', 'msavi_glcm_variance', 
                         'msavi_glcm_dissimilarity')

    if (!missing(dem)) {
        predictors <- stack(predictors, dem)
        predictor_names <- c(predictor_names, 'elev')
    }
    if (!missing(slopeaspect)) {
        predictors <- stack(predictors, slopeaspect$slope, aspect_cut)
        predictor_names <- c(predictor_names, 'slope', 'aspect')
    }

    predictors_filename <- file.path(output_path,
                                     paste0(image_basename, '_predictors.', 
                                            ext))

    names(predictors) <- predictor_names
    predictors <- mask(predictors, image_mask, maskvalue=1, 
                       filename=predictors_filename, 
                       overwrite=overwrite, datatype='INT2S')
    names(predictors) <- predictor_names

    # Save a copy of the original masks file along with the predictors file, so 
    # the masks can be easily located later.
    predictors_mask_filename <- file.path(output_path,
                                          paste0(image_basename, 
                                                 '_predictors_masks.', ext))
    mask_stack <- writeRaster(mask_stack, filename=predictors_mask_filename, 
                              overwrite=overwrite, 
                              datatype=dataType(mask_stack)[1])

    timer <- stop_timer(timer, label='Writing predictors')

    timer <- stop_timer(timer, label='Predictor calculation')

    return(predictors)
}
