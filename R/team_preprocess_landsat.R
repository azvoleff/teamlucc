#' Preprocess surface reflectance imagery from the Landsat CDR archive
#'
#' @export
#' @importFrom rgeos gContains gUnion
#' @param image_list a list of paths to Landsat CDR images that have been 
#' converted to bsq files with teampy.
#' @param dem a list of digital elevation models (DEMs) that (when mosaiced) 
#' covers the full extent of all the images in the image_list.
#' @param n_cpus the number of CPUs to use for processes that can run in 
#' parallel
#' @param cleartmp whether to clear temp files on each run through the loop
#' @examples
#' \dontrun{
#' image_list <- list('H:/Data/TEAM/VB/Rasters/Landsat/1986_037_LT5/proc/lndsr.LT50150531986037XXX18')
#' dem_list <- list('H:/Data/TEAM/VB/Rasters/DEM/ASTER/ASTGTM2_N09W084_dem.tif',
#'                  'H:/Data/TEAM/VB/Rasters/DEM/ASTER/ASTGTM2_N09W085_dem.tif',
#'                  'H:/Data/TEAM/VB/Rasters/DEM/ASTER/ASTGTM2_N10W084_dem.tif',
#'                  'H:/Data/TEAM/VB/Rasters/DEM/ASTER/ASTGTM2_N09W085_dem.tif')
#' team_preprocess(image_list, dem_list)
#' }
team_preprocess <- function(image_list, DEM, output_dir, n_cpus=2, 
                            cleartmp=FALSE) {
    ################################################################################
    # Verify extents and projections of images and DEMs match (DEM projection 
    # doesn't have to match image projection, but all DEMs must have the same 
    # projections, and all images must have the same projections).
    image_stacks <- lapply(image_list, stack)
    dem_rasts <- lapply(dem_list, raster)

    image_prj <- projection(image_stacks[[1]])
    if (any(lapply(image_stacks, projection) != image_prj)) {
        stop("each image in image_list must have the same projection")
    }
    dem_prj <- projection(dem_rasts[[1]])
    if (any(lapply(dem_rasts, projection) != dem_prj)) {
        stop("each DEM in dem_list must have the same projection")
    }

    #TODO: download ASTER that aligns with CDR SR imagery
    # plot(dem_mosaic_extent)
    # plot(image_extent_polys[[1]], add=TRUE)
    # plot(image_extent_polys[[2]], add=TRUE)
    # plot(image_extent_polys[[3]], add=TRUE)
    # Verify the combined extent of the DEMs in dem_list covers the full area 
    # of the images in image_list
    image_extent_polys <- lapply(image_stacks, get_extent_poly)
    # Make sure the DEM extents and image extent polys are in same projection 
    dem_extent_polys <- lapply(dem_rasts, function(rast) get_extent_poly(projectExtent(rast, CRS(image_prj))))
    dem_mosaic_extent <- dem_extent_polys[[1]]
    for (dem_extent_poly in dem_extent_polys[-1]) {
        dem_mosaic_extent <- gUnion(dem_mosaic_extent, dem_extent_poly)
    }
    extents_contained <- unlist(lapply(image_extent_polys, function(ext) gContains(dem_mosaic_extent, ext)))
    if (!any(extents_contained)) {
        warning("DEM does not fully cover extent of images in image_list")
    }

    ################################################################################
    # Mosaic DEMs
    print('Mosaicing DEMs...')
    trackTime(action='start')
    # See http://bit.ly/1dJPIeF re issue in raster that necessitates below workaround
    # TODO: Contact Hijmans re possible fix
    mosaicargs <- dem_rasts
    mosaicargs$fun <- mean
    dem_mosaic <- do.call(mosaic, mosaicargs)
    dem_mosaic_filename <- file.path(data_dir, 'Rasters/DEM/ASTER/ASTER_DEM_L5TSR_1986_mosaic.envi')
    dem_mosaic <- projectRaster(dem_mosaic, crs=CRS(image_prj), dem_mosaic_filename)
    #dem_mosaic <- mosaic(dem_rasts, fun='mean')
    trackTime()

    print('Running slopeasp_seq...')
    trackTime(action='start')
    slopeaspect_filename <- file.path(data_dir, 'Rasters/DEM/ASTER/ASTER_DEM_L5TSR_1986_slopeaspect.envi')
    slopeaspect <- slopeasp_seq(dem_mosaic, filename=slopeaspect_filename, overwrite=TRUE)
    trackTime()

    for (image_path in image_list) {
        ################################################################################
        # Load data and mask out clouds and missing values
        print('Loading data and masking clouds and missing data...')
        trackTime(action='start')
        image_rast_original_file <- file.path(data_dir, paste0(img_basename, '.bsq'))
        data_bands <- grep('band_[123457]_reflectance',  
                           get_band_names_from_hdr(extension(image_rast_original_file, '.hdr')))
        image_rast <- stack(image_rast_original_file, bands=data_bands)
        names(image_rast) <- get_band_names_from_hdr(extension(image_rast_original_file, '.hdr'))[data_bands]
        # Rename bands with shorter names. "sr" stands for "surface reflectance"
        names(image_rast) <- gsub('_reflectance', 'sr', names(image_rast))
        names(image_rast) <- gsub('band_', 'b', names(image_rast))

        mask_bands <- grep('(combined_cloud_mask)|(missing_mask)', 
                           get_band_names_from_hdr(extension(image_rast_original_file, '.hdr')))
        masks <- stack(image_rast_original_file, bands=mask_bands)
        names(masks) <- get_band_names_from_hdr(extension(image_rast_original_file, '.hdr'))[mask_bands]

        # The combined cloud mask includes the cloud_QA, cloud_shadow_QA, and 
        # adjacent_cloud_QA layers. Missing or clouded pixels are coded as 0, while 
        # good pixels are coded as 1.
        image_rast <- image_rast * masks$combined_cloud_mask * masks$missing_mask

        mask_file_path <- file.path(data_dir, paste0(img_basename, '_masked.envi'))
        writeRaster(image_rast, mask_file_path)
        trackTime()

        ################################################################################
        # Perform topographic correction
        print('Cropping slope/aspect raster to extent of Landsat image...')
        trackTime(action='start')
        beginCluster(2)
        matched_DEM_file <- file.path(data_dir, 'Rasters/DEM/ASTER/ASTER_DEM_L5TSR_1986.grd')
        dem_mosaic <- match_rasters(image_rast, dem_mosaic, filename=matched_DEM_file)
        endCluster()
        names(dem_mosaic) <- 'elevation'
        # Need to detach snow or else sfQuickInit will not work
        detach('package:snow', unload=TRUE)
        trackTime()

        print('Running topocorr...')
        library(spatial.tools)
        sfQuickInit(3)
        trackTime(action='start')
        image_rast_original_file <- file.path(data_dir, paste0(img_basename, '.bsq'))
        slopeaspect_filename <- file.path(data_dir, 'Rasters/DEM/ASTER/ASTER_DEM_L5TSR_1986_slopeaspect.grd')
        mask_file_path <- file.path(data_dir, paste0(img_basename, '_masked.envi'))
        slopeaspect <- brick(slopeaspect_filename)
        image_rast <- brick(mask_file_path)
        metadatafile <- extension(image_rast_original_file, 'txt')
        sunelev <- 90 - as.numeric(get_metadata_item(metadatafile, 'SolarZenith'))
        sunazimuth <- as.numeric(get_metadata_item(metadatafile, 'SolarAzimuth'))
        topocorr_filename <- file.path(data_dir, 'Rasters/Landsat/1986_037_LT5/proc/L5TSR_1986_topocorr.grd')
        # First draw a sample for the Minnaert k regression
        horizcells <- 10
        vertcells <- 10
        nsamp <- 200000 / (horizcells * vertcells)
        # Note that rowmajor indices are needed as raster layers are stored in 
        # rowmajor order, unlike most R objects that are addressed in column 
        # major order
        sampleindices <- gridsample(image_rast, horizcells=10, vertcells=10, 
                                    nsamp=nsamp, rowmajor=TRUE)
        image_rast <- topographic_corr(image_rast, slopeaspect, sunelev, 
                                       sunazimuth, method='minnaert_full', 
                                       filename=topocorr_filename, 
                                       inparallel=TRUE, overwrite=TRUE,
                                       sampleindices=sampleindices)
        trackTime()
        sfQuickStop()

        ################################################################################
        # Calculate additional predictor layers (MSAVI and textures)

        print('Calculating MSAVI2...')
        trackTime(action='start')
        MSAVI2_filename <- file.path(data_dir, 'Rasters/Landsat/1986_037_LT5/proc/L5TSR_1986_MSAVI.envi')
        MSAVI2_layer <- MSAVI2(red=raster(image_rast, layer=3),
                               nir=raster(image_rast, layer=4))
        writeRaster(MSAVI2_layer, MSAVI2_filename)
        trackTime()

        print('Calculating GLCM textures from MSAVI image...')
        trackTime(action='start')
        MSAVI2_glcm_filename <- file.path(data_dir, 'Rasters/Landsat/1986_037_LT5/proc/L5TSR_1986_MSAVI2_glcm.envi')
        MSAVI2_layer <- raster(MSAVI2_filename)
        MSAVI2_glcm <- glcm(MSAVI2_layer, filename=glcm_filename, overwrite=TRUE)
        trackTime()

        ################################################################################
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
                                  dem_mosaic,
                                  slopeaspect$slope,
                                  slopeaspect$aspect)
        image_rast_preds_filename <- file.path(data_dir, 'Rasters/Landsat/1986_037_LT5/proc/L5TSR_1986_preds.envi')
        writeRaster(image_rast_preds, image_rast_preds_filename)
        if (cleartmp) removeTmpFiles(h=1)
    }

}
