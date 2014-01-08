#' Preprocess surface reflectance imagery from the Landsat CDR archive
#'
#' @export
#' @import spatial.tools
#' @importFrom rgeos gContains gUnion
#' @param image_list a list of paths to a set of Landsat CDR images (for a 
#' single TEAM site) that have been converted to bsq files with teampy.
#' @param dem path to a digital elevation model (DEM) covering the full extent 
#' of all the images in \code{image_list}. See \code{team_setup_dem} for a 
#' function simplifying this.
#' @param slopeaspect path to a two layer raster stack with slope and aspect 
#' calculated from the above DEM
#' @param sitecode code to use as a prefix for all filenames
#' @param output_path the path to use for the output @param n_cpus the number 
#' of CPUs to use for processes that can run in parallel
#' @param cleartmp whether to clear temp files on each run through the loop
#' @param overwrite whether to overwrite existing files (otherwise an error 
#' will be raised)
#' @param notify notifier to use (defaults to \code{print} function). See the 
#' \code{notifyR} package for one way of sending notifications from R. The 
#' \code{notify} function should accept a string as the only argument.
#' @examples
#' \dontrun{
#' image_list <- c('H:/Data/TEAM/VB/Rasters/Landsat/1986_037_LT5/proc/lndsr.LT50150531986037XXX18.bsq',
#'                 'H:/Data/TEAM/VB/Rasters/Landsat/2001_014_LT5/proc/lndsr.LT50150532001014AAA01.bsq',
#'                 'H:/Data/TEAM/VB/Rasters/Landsat/2012_021_LE7/proc/lndsr.LE70150532012021EDC00.bsq')
#' dem <- 'H:/Data/TEAM/VB/LCLUC_Analysis/VB_dem_mosaic.envi'
#' slopeaspect <- 'H:/Data/TEAM/VB/LCLUC_Analysis/VB_dem_mosaic_slopeaspect.envi'
#' team_preprocess(image_list, dem, slopeaspect, "VB", 
#' 'H:/Data/TEAM/VB/LCLUC_Analysis', 3, TRUE)
#' }
team_preprocess_landsat <- function(image_list, dem, slopeaspect, sitecode, 
                                    output_path, n_cpus=1, cleartmp=FALSE, 
                                    overwrite=FALSE, notify=print) {
    notify("Starting preprocessing...")
    if (n_cpus > 1) sfQuickInit(n_cpus)
    ##########################################################################
    # Verify extents and projections of images and DEMs match (DEM projection 
    # doesn't have to match image projection, but all DEMs must have the same 
    # projections, and all images must have the same projections).
    image_stacks <- lapply(image_list, stack)

    image_prj <- projection(image_stacks[[1]])
    if (any(lapply(image_stacks, projection) != image_prj)) {
        stop("each image in image_list must have the same projection")
    }

    #TODO: download ASTER that aligns with CDR SR imagery
    # Verify the combined extent of the DEMs in dem_list covers the full area 
    # of the images in image_list
    image_extent_polys <- lapply(image_stacks, get_extent_poly)
    if ('RasterLayer' != class(dem)) dem <- raster(dem)
    # Make sure the DEM extents and image extent polys are in same projection 
    dem_extent_poly <- get_extent_poly(dem)
    extents_contained <- unlist(lapply(image_extent_polys,
                                       function(ext) gContains(dem_extent_poly, ext)))
    for (n in 1:length(extents_contained)) {
        warning(paste("DEM does not fully cover extent of", image_list[n]))
    }

    for (image_path in image_list) {
        notify(paste("Processing", image_path))

        ######################################################################
        # Determine image basename for use in naming subsequent files
        aq_date <- get_metadata_item(extension(image_path, '.txt'), 'AcquisitionDate')
        aq_date <- strptime(aq_date, format="%Y-%m-%dT%H:%M:%OSZ")
        short_name  <- get_metadata_item(extension(image_path, '.txt'), 'ShortName')
        image_basename <- paste(format(aq_date, '%Y_%j'), short_name, sep='_')

        ######################################################################
        # Load data and mask out clouds and missing values
        notify('Loading data and masking clouds and missing data...')
        notify(track_time(action='start'))
        data_bands <- grep('band_[123457]_reflectance',  
                           get_band_names_from_hdr(extension(image_path, '.hdr')))
        image_rast <- stack(image_path, bands=data_bands)
        names(image_rast) <- get_band_names_from_hdr(extension(image_path, '.hdr'))[data_bands]
        # Rename bands with shorter names. "sr" stands for "surface reflectance"
        names(image_rast) <- gsub('_reflectance', 'sr', names(image_rast))
        names(image_rast) <- gsub('band_', 'b', names(image_rast))

        mask_bands <- grep('(combined_cloud_mask)|(missing_mask)', 
                           get_band_names_from_hdr(extension(image_path, '.hdr')))
        masks <- stack(image_path, bands=mask_bands)
        names(masks) <- get_band_names_from_hdr(extension(image_path, '.hdr'))[mask_bands]

        # The combined cloud mask includes the cloud_QA, cloud_shadow_QA, and 
        # adjacent_cloud_QA layers. Missing or clouded pixels are coded as 0, while 
        # good pixels are coded as 1.
        image_rast_mask_path <- file.path(output_path,
                                          paste(sitecode, image_basename, 
                                                'mask.envi', sep='_'))
        image_rast_mask <- overlay(masks,
                                   fun=function(cloud_mask, missing_mask) {
                                       cloud_mask * missing_mask
                                   },
                                   filename=image_rast_mask_path, 
                                   overwrite=overwrite, 
                                   datatype=dataType(masks))

        image_rast_masked_path <- file.path(output_path,
                                            paste(sitecode, image_basename, 
                                                  'masked.envi', sep='_'))
        image_rast <- mask(image_rast, image_rast_mask, maskvalue=0, 
                           filename=image_rast_masked_path, 
                           overwrite=overwrite, datatype=dataType(image_rast))
        # image_rast_mask is no longer needed, so unload it to save memory
        rm(image_rast_mask)

        notify(track_time())

        ######################################################################
        # Perform topographic correction
        notify('Cropping dem and slope/aspect rasters to extent of Landsat image...')
        notify(track_time(action='start'))
        cropped_dem_file <- file.path(output_path,
                                      paste(sitecode, image_basename, 
                                            'dem.envi', sep='_'))
        cropped_dem <- match_rasters(image_rast, dem, 
                                     filename=cropped_dem_file, 
                                     overwrite=overwrite)
        if (!(class(slopeaspect) %in% c('RasterStack', 'RasterBrick'))) {
            slopeaspect <- stack(slopeaspect)
        }
        cropped_slopeaspect_file <- file.path(output_path,
                                              paste(sitecode, image_basename, 
                                                    'dem_slopeaspect.envi', 
                                                    sep='_'))
        #TODO: Check why extend in this match_rasters call is raising warning
        cropped_slopeaspect <- match_rasters(image_rast, slopeaspect, 
                                             filename=cropped_slopeaspect_file, 
                                             overwrite=overwrite)
        notify(track_time())
        if (cleartmp) removeTmpFiles(h=1)

        notify('Running topocorr...')
        notify(track_time(action='start'))

        # Draw a sample for the Minnaert k regression
        horizcells <- 10
        vertcells <- 10
        nsamp <- 200000 / (horizcells * vertcells)
        # Note that rowmajor indices are needed as raster layers are stored in 
        # rowmajor order, unlike most R objects that are addressed in column 
        # major order
        sampleindices <- gridsample(image_rast, horizcells=10, vertcells=10, 
                                    nsamp=nsamp, rowmajor=TRUE)

        metadatafile <- extension(image_path, 'txt')
        sunelev <- 90 - as.numeric(get_metadata_item(metadatafile, 'SolarZenith'))
        sunazimuth <- as.numeric(get_metadata_item(metadatafile, 'SolarAzimuth'))

        topocorr_filename <- file.path(output_path,
                                       paste(sitecode, image_basename, 
                                             'masked_tc.envi', sep='_'))
        if (n_cpus > 1) {
            inparallel <- TRUE
        } else {
            inparallel <- FALSE
        }
        image_rast <- topographic_corr(image_rast, cropped_slopeaspect, sunelev, 
                                       sunazimuth, method='minnaert_full', 
                                       filename=topocorr_filename, 
                                       inparallel=inparallel, 
                                       overwrite=overwrite, 
                                       sampleindices=sampleindices)
        notify(track_time())
        if (cleartmp) removeTmpFiles(h=1)

        ######################################################################
        # Calculate additional predictor layers (MSAVI and textures)
        notify('Calculating MSAVI2...')
        notify(track_time(action='start'))
        MSAVI2_filename <- file.path(output_path,
                                     paste(sitecode, image_basename, 
                                           'masked_tc_MSAVI2.envi', sep='_'))
        MSAVI2_layer <- MSAVI2(red=raster(image_rast, layer=3),
                               nir=raster(image_rast, layer=4))
        MSAVI2_layer <- writeRaster(MSAVI2_layer, MSAVI2_filename, 
                                    overwrite=overwrite, datatype=dataType(MSAVI2_layer))
        notify(track_time())

        notify('Calculating GLCM textures from MSAVI image...')
        notify(track_time(action='start'))
        MSAVI2_glcm_filename <- file.path(output_path,
                                          paste(sitecode, image_basename, 
                                                'masked_tc_MSAVI2_glcm.envi', 
                                                sep='_'))
        MSAVI2_layer <- raster(MSAVI2_filename)
        MSAVI2_glcm <- glcm(MSAVI2_layer)
        MSAVI2_glcm <- writeRaster(MSAVI2_glcm, filename=MSAVI2_glcm_filename, 
                    overwrite=overwrite, datatype=dataType(MSAVI2_glcm))
        notify(track_time())

        ######################################################################
        # Layer stack predictor layers:
        image_rast_preds <- stack(raster(image_rast, layer=1),
                                  raster(image_rast, layer=2),
                                  raster(image_rast, layer=3),
                                  raster(image_rast, layer=4),
                                  raster(image_rast, layer=5),
                                  raster(image_rast, layer=6),
                                  MSAVI2_layer,
                                  MSAVI2_glcm$glcm_mean,
                                  MSAVI2_glcm$glcm_variance,
                                  MSAVI2_glcm$glcm_dissimilarity,
                                  cropped_dem,
                                  cropped_slopeaspect$slope,
                                  cropped_slopeaspect$aspect)
        image_rast_preds_filename <- file.path(output_path,
                                               paste(sitecode, image_basename, 
                                                     'predictors.envi', 
                                                     sep='_'))
        image_rast_preds <- writeRaster(image_rast_preds, 
                                        image_rast_preds_filename, 
                                        overwrite=overwrite, 
                                        datatype=dataType(image_rast_preds))
        if (cleartmp) removeTmpFiles(h=1)
    }
    if (n_cpus > 1) sfQuickStop()
}
