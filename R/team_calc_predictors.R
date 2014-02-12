#' Calculate predictor layers for a classification
#'
#' @export
#' @importFrom glcm glcm
#' @param image_dirs list of paths to a set of Landsat images that have been 
#' preprocessed by the \code{team_preprocess_landsat} function.
#' @param dem_path path to a set of DEMs as output by \code{team_setup_dem}
#' @param output_path the path to use for the output
#' @param n_cpus the number of CPUs to use for processes that can run in 
#' parallel
#' @param cleartmp whether to clear temp files on each run through the loop
#' @param overwrite whether to overwrite existing files (otherwise an error 
#' will be raised)
#' @param notify notifier to use (defaults to \code{print} function). See the 
#' \code{notifyR} package for one way of sending notifications from R. The 
#' \code{notify} function should accept a string as the only argument.
team_calc_predictors <- function(image_dirs, dem_path, output_path, n_cpus=1, 
                                 cleartmp=FALSE, overwrite=FALSE,
                                 notify=print) {
    if (!file_test("-d", dem_path)) {
        stop(paste(dem_path, "does not exist"))
    }
    if (!file_test("-d", output_path)) {
        stop(paste(output_path, "does not exist"))
    }

    timer <- Track_time(notify)
    timer <- start_timer(timer, label='Predictor calculation')

    if (n_cpus > 1) beginCluster(n_cpus)

    # Setup a regex to identify preprocessed images
    preproc_regex <- '^{a-zA-Z}{2,3}_[0-9]{6}_[0-9]_L[457][ET]SR_masked_tc'

    for (image_dir in image_dirs) {
        if (!file_test("-d", image_dir)) {
            stop(paste(image_dir, "does not exist"))
        }
        preproc_files <- dir(image_dir, pattern=preproc_regex)
        image_basenames <- unique(str_extract(preproc_files, preproc_regex))
        if (length(image_basenames) == 0) {
            stop(paste('no files found in', image_dir))
        }

        for (image_basename in image_basenames) {
            # Choose between cloud filled, gap filled, or original files, 
            # depending on what is available
            poss_files <- c(file.path(image_dir, paste0(image_basename, 
                                                        '_cf_gf.envi')),
                            file.path(image_dir, paste0(image_basename, 
                                                        '_gf.envi')),
                            file.path(image_dir, paste0(image_basename, 
                                                        '_cf.envi')),
                            file.path(image_dir, paste0(image_basename, 
                                                        '.envi')))
            exist_files <- unlist(lapply(poss_files, file.exists))
            if (sum(exist_files) > 1) {
                stop(paste('multiple file versions exist for', image_basename))
            } else if (sum(exist_files) == 0) {
                stop(paste('no pre-processed files found for', image_basename))
            }

            image_path <- poss_files[exist_files]

            # Update the basename to refer the chosen file
            image_basename <- basename(file_path_sans_ext(image_path))

            image_stack <- brick(image_path)

            ######################################################################
            # Calculate additional predictor layers (MSAVI and textures)
            timer <- start_timer(timer, label=paste(image_basename, '-', 'MSAVI2'))
            MSAVI2_filename <- file.path(output_path,
                                         paste(image_basename, 'MSAVI2.envi', 
                                               sep='_'))
            MSAVI2_layer <- MSAVI2(red=raster(image_stack, layer=3),
                                   nir=raster(image_stack, layer=4))
            # Truncate MSAVI2 to range between 0 and 1, and scale by 10,000 so it 
            # can be saved as a INT2S
            MSAVI2_layer <- calc(MSAVI2_layer, fun=function(vals) {
                    vals[vals > 1] <- 1
                    vals[vals < 0] <- 0
                    vals <- round(vals * 10000)
                }, filename=MSAVI2_filename, overwrite=overwrite, datatype="INT2S")
            timer <- stop_timer(timer, label=paste(image_basename, '-', 'MSAVI2'))

            timer <- start_timer(timer, label=paste(image_basename, '-', 'glcm'))
            MSAVI2_glcm_filename <- file.path(output_path,
                                              paste(image_basename, 
                                                    'MSAVI2_glcm.envi', 
                                                    sep='_'))
            glcm_statistics <- c('mean', 'variance', 'homogeneity', 'contrast', 
                                 'dissimilarity', 'entropy', 'second_moment', 
                                 'correlation')
            # Note the min_x and max_x are given for MSAVI2 that has been scaled by 
            # 10,000
            MSAVI2_glcm <- apply_windowed(MSAVI2_layer, glcm, edge=c(1, 3), 
                                          min_x=0, max_x=10000, 
                                          filename=MSAVI2_glcm_filename, 
                                          overwrite=overwrite, 
                                          statistics=glcm_statistics, na_opt='center')
            names(MSAVI2_glcm) <- paste('glcm', glcm_statistics, sep='_')
            timer <- stop_timer(timer, label=paste(image_basename, '-', 'glcm'))

            ######################################################################
            # Load DEM, slope, and aspect, and reclass aspect
            timer <- start_timer(timer, label=paste(image_basename, '-', 'process dem and slopeaspect'))
            
            #TODO: load dem and slopeaspect layers
            cropped_dem <- NULL
            slopeaspect <- NULL
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
            timer <- stop_timer(timer, label=paste(image_basename, '-', 'process dem and slopeaspect'))

            ######################################################################
            # Layer stack predictor layers:
            timer <- start_timer(timer, label=paste(image_basename, '-', 'write predictors'))
            predictors <- stack(raster(image_stack, layer=1),
                                raster(image_stack, layer=2),
                                raster(image_stack, layer=3),
                                raster(image_stack, layer=4),
                                raster(image_stack, layer=5),
                                raster(image_stack, layer=6),
                                MSAVI2_layer,
                                scale_raster(MSAVI2_glcm$glcm_mean),
                                scale_raster(MSAVI2_glcm$glcm_variance),
                                scale_raster(MSAVI2_glcm$glcm_dissimilarity),
                                cropped_dem,
                                slopeaspect$slope,
                                aspect_cut)
            predictors_filename <- file.path(output_path,
                                             paste(image_basename, 
                                                   'predictors.envi', sep='_'))
            #TODO: load mask
            image_stack_mask <- NULL

            predictors <- mask(predictors, image_stack_mask, maskvalue=0, 
                               filename=predictors_filename, overwrite=overwrite,
                               datatype='INT2S')
            names(predictors) <- c('b1', 'b2', 'b3', 'b4', 'b5', 'b7', 'msavi', 
                                  'msavi_glcm_mean', 'msavi_glcm_variance', 
                                  'msavi_glcm_dissimilarity', 'elev', 'slope', 
                                  'aspect')
            timer <- stop_timer(timer, label=paste(image_basename, '-', 'write predictors'))

            timer <- stop_timer(timer, label=paste('Calculating predictors', image_basename))

            if (cleartmp) removeTmpFiles(h=1)
        }
    }
    if (n_cpus > 1) endCluster()

    timer <- stop_timer(timer, label='Predictor calculation')
}
