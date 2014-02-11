team_preprocess_landsat <- function(image_dirs, output_path, n_cpus=1, 
                                    cleartmp=FALSE, overwrite=FALSE, 
                                    notify=print) {

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

            image_stack <- brick(poss_files[exist_files])

            sitecode <- str_extract(image_basename, '^[a-zA-Z]{2,3}')

            ######################################################################
            # Calculate additional predictor layers (MSAVI and textures)
            timer <- start_timer(timer, label=paste(image_basename, '-', 'MSAVI2'))
            MSAVI2_filename <- file.path(output_path,
                                         paste(sitecode, image_basename, 
                                               'masked_tc_MSAVI2.envi', sep='_'))
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
                                              paste(sitecode, image_basename, 
                                                    'masked_tc_MSAVI2_glcm.envi', 
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
                                             paste(sitecode, image_basename, 
                                                   'predictors.envi', sep='_'))
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
