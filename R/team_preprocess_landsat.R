#' Preprocess surface reflectance imagery from the Landsat CDR archive
#'
#' This function preprocesses surface reflectance imagery from the Landsat 
#' Climate Data Record (CDR) archive. \code{team_preprocess_landsat} can 
#' reproject CDR tiles to match the projection of a given \code{aoi}, crop the 
#' tiles to match the \code{aoi} or a common WRS-2 path/row polygon, mask 
#' missing data and clouds out of the CDR tiles, and perform topographic 
#' correction.
#'
#' Prior to running \code{team_preprocess_landsat}, the 
#' \code{\link{espa_extract}} and should be used to extract the original 
#' zipfiles supplied by USGS.  \code{\link{unstack_ledapscdr}} should then be 
#' used to unstack the HDF format files and convert them into ENVI binary 
#' format. To perform topographic correction with 
#' \code{team_preprocess_landsat}, first run \code{\link{team_setup_dem}} to 
#' preprocess a set of DEM tiles. Then run \code{team_preprocess_landsat} with 
#' the \code{tc=TRUE} option.
#' @export
#' @importFrom rgeos gIntersection
#' @importFrom wrspathrow pathrow_poly
#' @importFrom tools file_path_sans_ext
#' @param image_dirs list of paths to a set of Landsat CDR image files in ENVI 
#' format as output by the \code{unstack_ledapscdr} function.
#' @param output_path the path to use for the output (optional - if NULL then 
#' output images will be saved alongside the input images in the same folder).
#' @param prefix string to use as a prefix for all filenames
#' @param tc whether to topographically correct imagery
#' @param aoi area of interest (AOI), as a \code{SpatialPolygonsDataFrame}, to 
#' use as as bounding box when selecting DEMs
#' @param dem_path path to a set of DEMs as output by \code{team_setup_dem} 
#' (only required if tc=TRUE)
#' @param n_cpus the number of CPUs to use for processes that can run in 
#' parallel
#' @param cleartmp whether to clear temp files on each run through the loop
#' @param overwrite whether to overwrite existing files (otherwise an error 
#' will be raised)
#' @param notify notifier to use (defaults to \code{print} function). See the 
#' \code{notifyR} package for one way of sending notifications from R. The 
#' \code{notify} function should accept a string as the only argument.
#' @param verbose whether to print detailed status messages and timing 
#' information
#' @seealso \code{\link{espa_extract}}, \code{\link{unstack_ledapscdr}}, 
#' \code{\link{team_setup_dem}}
team_preprocess_landsat <- function(image_dirs, prefix, tc=FALSE, aoi=NULL,
                                    dem_path=NULL, output_path=NULL, n_cpus=1, 
                                    cleartmp=FALSE,  overwrite=FALSE, 
                                    notify=print, verbose=FALSE) {
    if (tc && !file_test("-d", dem_path)) {
        stop(paste(dem_path, "does not exist"))
    }
    if (!is.null(output_path) && !file_test("-d", output_path)) {
        stop(paste(output_path, "does not exist"))
    }

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
            gain(image_stack) <- 1/as.numeric(get_metadata_item(band_files[1], 
                                                                'scale_factor'))
            offs(image_stack) <- as.numeric(get_metadata_item(band_files[1], 
                                                              'add_offset'))

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

    if (n_cpus > 1) endCluster()
    for (n in 1:length(image_stacks)) {
        if (n_cpus > 1) beginCluster(n_cpus)
        image_stack <- image_stacks[[n]]
        mask_stack <- mask_stacks[[n]]
        band1_imagefile <- image_files[[n]][1]

        if (is.null(output_path)) {
            this_output_path <- dirname(band1_imagefile)
        } else {
            this_output_path  <- output_path
        }

        ######################################################################
        # Determine image basename for use in naming subsequent files
        aq_date <- get_metadata_item(band1_imagefile, 'AcquisitionDate')
        aq_date <- strptime(aq_date, format="%Y-%m-%dT%H:%M:%OSZ")
        short_name  <- get_metadata_item(band1_imagefile, 'ShortName')
        WRS_Path <- sprintf('%03i', as.numeric(get_metadata_item(band1_imagefile, 'WRS_Path')))
        WRS_Row <- sprintf('%03i', as.numeric(get_metadata_item(band1_imagefile, 'WRS_Row')))
        image_basename <- paste0(WRS_Path, '-', WRS_Row, '_',
                                 format(aq_date, '%Y-%j'), '_', short_name)

        timer <- start_timer(timer, label=paste('Preprocessing', image_basename))

        #######################################################################
        # Reproject images to match the projection being used for this image.  
        # This is either the projection of aoi (if aoi is supplied), or the UTM 
        # zone of the centroid of this path and row.
        this_pathrow_poly <- pathrow_poly(as.numeric(WRS_Path), 
                                          as.numeric(WRS_Row))
        if (!is.null(aoi)) {
            to_proj4string <- proj4string(aoi)
        } else {
            to_proj4string <- utm_zone(this_pathrow_poly, proj4string=TRUE)
        }
        if (!proj4comp(to_proj4string, proj4string(image_stack))) {
            if (verbose) timer <- start_timer(timer,
                                              label=paste(image_basename, '- reproject'))
            image_stack <- projectRaster(image_stack, crs=CRS(to_proj4string))
            mask_stack <- projectRaster(mask_stack, crs=CRS(to_proj4string), method='ngb')
            if (verbose) timer <- stop_timer(timer,
                                             label=paste(image_basename, '- reproject'))
        } else {
            proj4string(image_stack) <- to_proj4string
            proj4string(mask_stack) <- to_proj4string
        }

        ######################################################################
        # Crop image as necessary
        if (verbose) timer <- start_timer(timer, label=paste(image_basename, 
                                                            '-', 'cropping'))
        # Note that the image_stack may have been reprojected to match the 
        # projection system of the AOI. Ensure the pathrow_poly projection 
        # system still matches:
        this_pathrow_poly <- spTransform(this_pathrow_poly, 
                                         CRS(proj4string(image_stack)))

        if (!is.null(aoi)) {
            # If an aoi IS supplied, match the image extent to that of the AOI 
            # cropped to the appropriate Landsat path/row polygon.
            crop_area <- gIntersection(this_pathrow_poly, aoi, byid=TRUE)
        } else {
            # If an aoi IS NOT supplied, match the image extent to the 
            # appropriate Landsat path/row polygon.
            crop_area <- this_pathrow_poly
        }
        image_stack <- crop(image_stack, crop_area)
        image_stack <- extend(image_stack, crop_area)
        mask_stack <- crop(mask_stack, crop_area)
        mask_stack <- extend(mask_stack, crop_area)
        if (verbose) timer <- stop_timer(timer, label=paste(image_basename, 
                                                            '-', 'cropping'))

        ######################################################################
        # Mask out clouds and missing values
        if (verbose) timer <- start_timer(timer, label=paste(image_basename, 
                                                             '-', 'masking'))

        # The combined cloud mask includes the cloud_QA, cloud_shadow_QA, and 
        # adjacent_cloud_QA layers. Missing or clouded pixels are coded as 0, while 
        # good pixels are coded as 1.
        image_stack_mask_path <- file.path(this_output_path,
                                           paste(prefix, image_basename, 
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
        image_stack <- mask(image_stack, image_stack_mask, maskvalue=0)

        mask_stack_path <- file.path(this_output_path,
                                     paste(prefix, image_basename, 
                                           'masks.envi', sep='_'))
        mask_stack <- writeRaster(mask_stack, filename=mask_stack_path, 
                                  overwrite=overwrite, 
                                  datatype='INT2S')
        if (verbose) timer <- stop_timer(timer, label=paste(image_basename, 
                                                            '-', 'masking'))

        ######################################################################
        # Perform topographic correction if tc=TRUE
        if (tc) {
            if (verbose) timer <- start_timer(timer, label=paste(image_basename, 
                                                                 '-', 'topocorr'))

            ######################################################################
            # Load dem, slope, and aspect
            slopeaspect_filename <- file.path(dem_path,
                                              paste0('slopeaspect_', 
                                                     WRS_Path, '-', WRS_Row, '.envi'))
            slopeaspect <- brick(slopeaspect_filename)

            if (!proj4comp(proj4string(image_stack), proj4string(slopeaspect))) {
                stop(paste0('slopeaspect and image_stack projections do not match.\nslopeaspect proj4string: ', 
                            proj4string(slopeaspect), '\nimage_stack proj4string: ',
                            proj4string(image_stack)))
            } else {
                # Projection strings are functionally identical - so make sure 
                # their textual representations are the same.
                proj4string(slopeaspect) <- proj4string(image_stack)
            }

            if (ncell(image_stack) > 400000) {
                # Draw a sample for the Minnaert k regression
                horizcells <- 10
                vertcells <- 10
                nsamp <- 200000 / (horizcells * vertcells)
                # Note that rowmajor indices are needed as raster layers are stored in 
                # rowmajor order, unlike most R objects that are addressed in column 
                # major order
                sampleindices <- gridsample(image_stack, horizcells=10, vertcells=10, 
                                            nsamp=nsamp, rowmajor=TRUE)
            } else {
                sampleindices <- NULL
            }
            # Remember that slopeaspect layers are scaled to INT2S, but 
            # topographic_corr expects them as floats, so apply the scale factors 
            # used in team_setup_dem
            slopeaspect_flt <- stack(raster(slopeaspect, layer=1) / 10000,
                                     raster(slopeaspect, layer=2) / 1000)
            sunelev <- 90 - as.numeric(get_metadata_item(band1_imagefile, 'SolarZenith'))
            sunazimuth <- as.numeric(get_metadata_item(band1_imagefile, 'SolarAzimuth'))
            output_filename <- file.path(this_output_path,
                                         paste(prefix, image_basename, 
                                               'tc.envi', sep='_'))
            image_stack <- topographic_corr(image_stack, slopeaspect_flt, sunelev, 
                                            sunazimuth, method='minnaert_full', 
                                            asinteger=TRUE, 
                                            sampleindices=sampleindices)
            if (verbose) timer <- stop_timer(timer, label=paste(image_basename, 
                                                                '-', 'topocorr'))

        } else {
            # Skip topographic correction, so don't append _tc to filename
            output_filename <- file.path(this_output_path,
                                         paste0(prefix, '_', image_basename, 
                                                '.envi'))
        }

        image_stack <- writeRaster(image_stack, filename=output_filename, 
                                   overwrite=overwrite, datatype='INT2S')

        timer <- stop_timer(timer, label=paste('Preprocessing', image_basename))

        if (cleartmp) removeTmpFiles(h=1)
        if (n_cpus > 1) endCluster()
    }

    timer <- stop_timer(timer, label='Preprocessing images')
}
