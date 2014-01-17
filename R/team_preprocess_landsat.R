#' Preprocess surface reflectance imagery from the Landsat CDR archive
#'
#' @export
#' @importFrom rgeos gContains gUnion
#' @param image_dirs list of paths to a set of Landsat CDR image files in ENVI 
#' format as output by the \code{unstack_ledapscdr} function.
#' @param dem path to a digital elevation model (DEM) covering the full extent 
#' of all the images in \code{image_dirs}. See \code{team_setup_dem} for a 
#' function simplifying this.
#' @param slopeaspect path to a two band raster of slope (band 1) and aspect 
#' (band 2) in a format readable by \code{raster}
#' @param sitecode code to use as a prefix for all filenames
#' @param output_path the path to use for the output
#' @param aoi an area of interest (AOI) to crop from each image
#' @param n_cpus the number of CPUs to use for processes that can run in 
#' parallel
#' @param cleartmp whether to clear temp files on each run through the loop
#' @param overwrite whether to overwrite existing files (otherwise an error 
#' will be raised)
#' @param notify notifier to use (defaults to \code{print} function). See the 
#' \code{notifyR} package for one way of sending notifications from R. The 
#' \code{notify} function should accept a string as the only argument.
#' @examples
#' \dontrun{
#' image_dirs <- c('H:/Data/TEAM/VB/Rasters/Landsat/1986_037_LT5/proc',
#'                 'H:/Data/TEAM/VB/Rasters/Landsat/2001_014_LT5/proc',
#'                 'H:/Data/TEAM/VB/Rasters/Landsat/2012_021_LE7/proc')
#' dem <- 'H:/Data/TEAM/VB/LCLUC_Analysis/VB_dem_mosaic.envi'
#' slopeaspect <- 'H:/Data/TEAM/VB/LCLUC_Analysis/VB_dem_mosaic_slopeaspect.envi'
#' team_preprocess(image_dirs, dem, slopeaspect, "VB", 
#'                 'H:/Data/TEAM/VB/LCLUC_Analysis', 3, TRUE)
#' }
team_preprocess_landsat <- function(image_dirs, dem, slopeaspect, sitecode, 
                                    output_path, aoi=NULL, n_cpus=1, 
                                    cleartmp=FALSE,  overwrite=FALSE, notify=print) {
    timer <- Track_time(notify)

    timer <- start_timer(timer, label='Preprocessing images')
    if (n_cpus > 1) beginCluster(n_cpus)

    # Setup a regex to identify Landsat CDR images
    lndsr_regex <- '^lndsr.((LT4)|(LT5)|(LE7)|(LE8))[0-9]{6}[12][0-9]{6}[a-zA-Z]{3}[0-9]{2}'

    image_bands <- c('band1', 'band2', 'band3', 'band4', 'band5', 'band7')
    mask_bands <- c('fill_QA', 'fmask_band')
    image_files <- c()
    image_stacks <- c()
    mask_stacks <- c()
    for (image_dir in image_dirs) {
        if (!file_test("-d", image_dir)) {
            stop(paste(image_dir, "does not exist"))
        }
        lndsr_files <- dir(image_dir, pattern=lndsr_regex)
        loc <- regexpr(lndsr_regex, lndsr_files)
        stop_char <- loc + attr(loc, 'match.length') - 1
        image_basenames <- unique(substr(lndsr_files, loc, stop_char))

        if (length(image_basenames) == 0) {
            stop(paste('no files found in', image_dir))
        }

        for (image_basename in image_basenames) {
            band_files <- c()
            for (image_band in image_bands) {
                band_files <- c(band_files,
                                paste(file.path(image_dir, image_basename), 
                                      image_band, sep='_'))
            }
            band_files <- paste0(band_files, '.envi')
            image_stack <- stack(band_files)
            names(image_stack) <- image_bands
            image_files <- c(image_files, list(band_files))
            # Read gain and offset from the metadata file as raster does not 
            # support passing as.is=TRUE to readGDAL, meaning that GDAL 
            # automatically scales the returned data turning it into floats.
            gain(image_stack) <- 1/as.numeric(get_metadata_item(band_files[1], 'scale_factor'))
            offs(image_stack) <- as.numeric(get_metadata_item(band_files[1], 'add_offset'))

            mask_band_files <- c()
            for (mask_band in mask_bands) {
                mask_band_files <- c(mask_band_files,
                                     paste(file.path(image_dir, 
                                                     image_basename), 
                                           mask_band, sep='_'))
            }
            mask_band_files <- paste0(mask_band_files, '.envi')
            mask_stack <- stack(mask_band_files)
            names(mask_stack) <- mask_bands

            image_stacks <- c(image_stacks, image_stack)
            mask_stacks <- c(mask_stacks, mask_stack)
        }
    }

    ##########################################################################
    # Verify extents and projections of images and DEMs match (DEM projection 
    # doesn't have to match image projection, but all DEMs must have the same 
    # projections, and all images must have the same projections).
    image_prj <- projection(image_stacks[[1]])
    if (any(lapply(image_stacks, projection) != image_prj)) {
        stop("each input image must have the same projection")
    }

    # Verify the combined extent of the DEMs in dem_list covers the full area 
    # of the images in image_stacks
    image_extent_polys <- lapply(image_stacks, get_extent_poly)
    if ('RasterLayer' != class(dem)) dem <- raster(dem)
    # Make sure the DEM extents and image extent polys are in same projection 
    dem_extent_poly <- get_extent_poly(dem)
    extents_contained <- unlist(lapply(image_extent_polys,
                                       function(ext) gContains(dem_extent_poly, ext)))
    for (n in 1:length(extents_contained)) {
        warning(paste("DEM does not fully cover extent of", image_files[n]))
    }

    for (n in 1:length(image_stacks)) {
        image_stack <- image_stacks[[n]]
        mask_stack <- mask_stacks[[n]]
        band1_imagefile <- image_files[[n]][1]

        ######################################################################
        # Determine image basename for use in naming subsequent files
        aq_date <- get_metadata_item(band1_imagefile, 'AcquisitionDate')
        aq_date <- strptime(aq_date, format="%Y-%m-%dT%H:%M:%OSZ")
        short_name  <- get_metadata_item(band1_imagefile, 'ShortName')
        WRS_Path <- sprintf('%03i', as.numeric(get_metadata_item(band1_imagefile, 'WRS_Path')))
        WRS_Row <- sprintf('%03i', as.numeric(get_metadata_item(band1_imagefile, 'WRS_Row')))
        image_basename <- paste(paste0(WRS_Path, WRS_Row),
                                format(aq_date, '%Y%j'), short_name, sep='_')
        if (!is.null(aoi)) {
            image_basename <- paste(image_basename, 'crop', sep='_')
        }

        timer <- start_timer(timer, label=paste('Preprocessing', image_basename))

        ######################################################################
        # Crop image to AOI if desired
        if (!is.null(aoi)) {
            timer <- start_timer(timer, label=paste(image_basename, '-', 'crop'))
            if (class(aoi) != 'SpatialPolygonsDataFrame') {
                stop('aoi must be a SpatialPolygonsDataFrame')
            } else if (projection(aoi) != projection(image_stack)) {
                stop(paste('projections of aoi and', image_basename, 'do not match'))
            }
            image_stack <- crop(image_stack, aoi)
            mask_stack <- crop(mask_stack, aoi)
            timer <- stop_timer(timer, label=paste(image_basename, '-', 'crop'))
        }

        ######################################################################
        # Load data and mask out clouds and missing values
        timer <- start_timer(timer, label=paste(image_basename, '-', 'masking'))

        # The combined cloud mask includes the cloud_QA, cloud_shadow_QA, and 
        # adjacent_cloud_QA layers. Missing or clouded pixels are coded as 0, while 
        # good pixels are coded as 1.
        image_stack_mask_path <- file.path(output_path,
                                          paste(sitecode, image_basename, 
                                                'mask.envi', sep='_'))
        # fmask_band key:
        # 	0 = clear
        # 	1 = water
        # 	2 = cloud_shadow
        # 	3 = snow
        # 	4 = cloud
        # 	255 = fill value
        # fill_QA key:
        # 	0 = not fill
        # 	255 = fill
        image_stack_mask <- overlay(mask_stack$fmask_band, mask_stack$fill_QA,
            fun=function(fmask, fill) {
                ((fmask == 0) | (fmask == 1)) & (fill == 0)
                })

        image_stack_masked_path <- file.path(output_path,
                                            paste(sitecode, image_basename, 
                                                  'masked.envi', sep='_'))
        image_stack <- mask(image_stack, image_stack_mask, maskvalue=0,
                            filename=image_stack_masked_path, 
                            overwrite=overwrite, datatype=dataType(image_stack)[1])
        timer <- stop_timer(timer, label=paste(image_basename, '-', 'masking'))

        ######################################################################
        # Crop dem, slope, and aspect
        timer <- start_timer(timer, label=paste(image_basename, '-', 'crop DEM'))
        cropped_dem_file <- file.path(output_path,
                                      paste(sitecode, image_basename, 
                                            'dem.envi', sep='_'))
        cropped_dem <- match_rasters(image_stack, dem)
        cropped_dem <- calc(cropped_dem, fun=function(vals) {
                round(vals)
            }, filename=cropped_dem_file, overwrite=overwrite, 
            datatype='INT2S')
        names(cropped_dem) <- "dem"
        timer <- stop_timer(timer, label=paste(image_basename, '-', 'crop DEM'))

        timer <- start_timer(timer, label=paste(image_basename, '-', 'calculate slope/aspect'))
        slopeaspect_cropped_file <- file.path(output_path,
                                      paste(sitecode, image_basename, 
                                            'dem_slopeaspect.envi', sep='_'))
        if (!(class(slopeaspect) %in% c('RasterStack', 'RasterBrick'))) {
            slopeaspect <- brick(slopeaspect)
        }
        slopeaspect <- match_rasters(image_stack, slopeaspect, 
                                     filename=slopeaspect_cropped_file, 
                                     overwrite=overwrite)
        names(slopeaspect) <- c('slope', 'aspect')
        timer <- stop_timer(timer, label=paste(image_basename, '-', 'calculate slope/aspect'))

        ######################################################################
        # Perform topographic correction
        timer <- start_timer(timer, label=paste(image_basename, '-', 'topocorr'))
        # Draw a sample for the Minnaert k regression
        horizcells <- 10
        vertcells <- 10
        nsamp <- 200000 / (horizcells * vertcells)
        # Note that rowmajor indices are needed as raster layers are stored in 
        # rowmajor order, unlike most R objects that are addressed in column 
        # major order
        sampleindices <- gridsample(image_stack, horizcells=10, vertcells=10, 
                                    nsamp=nsamp, rowmajor=TRUE)

        sunelev <- 90 - as.numeric(get_metadata_item(band1_imagefile, 'SolarZenith'))
        sunazimuth <- as.numeric(get_metadata_item(band1_imagefile, 'SolarAzimuth'))

        topocorr_filename <- file.path(output_path,
                                       paste(sitecode, image_basename, 
                                             'masked_tc.envi', sep='_'))
        if (n_cpus > 1) {
            inparallel <- TRUE
        } else {
            inparallel <- FALSE
        }
        image_stack <- topographic_corr(image_stack, slopeaspect, sunelev, 
                                       sunazimuth, method='minnaert_full', 
                                       filename=topocorr_filename, 
                                       inparallel=inparallel, 
                                       overwrite=overwrite, datatype='INT2S',
                                       sampleindices=sampleindices,
                                       scale_factor=10000, asinteger=TRUE)
        timer <- stop_timer(timer, label=paste(image_basename, '-', 'topocorr'))

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
                                      statistics=glcm_statistics)
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
                            MSAVI2_glcm$glcm_mean,
                            MSAVI2_glcm$glcm_variance,
                            MSAVI2_glcm$glcm_dissimilarity,
                            cropped_dem$dem,
                            slopeaspect$slope,
                            slopeaspect$aspect)
        predictors_filename <- file.path(output_path,
                                         paste(sitecode, image_basename, 
                                               'predictors.envi', sep='_'))
        predictors <- writeRaster(predictors, filename=predictors_filename, 
                                  overwrite=overwrite)
        names(predictors) <- c('b1', 'b2', 'b3', 'b4', 'b5', 'b7', 'msavi', 
                              'msavi_glcm_mean', 'msavi_glcm_variance', 
                              'msavi_glcm_dissimilarity', 'elev', 'slope', 
                              'aspect')
        timer <- stop_timer(timer, label=paste(image_basename, '-', 'write predictors'))

        timer <- stop_timer(timer, label=paste('Preprocessing', image_basename))

        if (cleartmp) removeTmpFiles(h=1)
    }
    if (n_cpus > 1) endCluster()

    timer <- stop_timer(timer, label='Preprocessing images')
}
