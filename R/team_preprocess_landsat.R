#' Preprocess surface reflectance imagery from the Landsat CDR archive
#'
#' @export
#' @importFrom wrspathrow pathrow_poly
#' @importFrom rgdal readOGR
#' @importFrom rgeos gContains gIntersection
#' @importFrom tools file_path_sans_ext
#' @param image_dirs list of paths to a set of Landsat CDR image files in ENVI 
#' format as output by the \code{unstack_ledapscdr} function.
#' @param output_path the path to use for the output
#' @param dem_path path to a set of DEMs as output by \code{team_setup_dem}
#' @param sitecode code to use as a prefix for all filenames
#' @param aoi_file an area of interest (AOI) to crop from each image, in a file 
#' format readable by \code{readOGR}
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
#' team_preprocess(image_dirs, 'H:/Data/TEAM/VB/LCLUC_Analysis', 'VB'
#' 'H:/Data/TEAM/VB/LCLUC_Analysis', n_cpus=3)
#' }
team_preprocess_landsat <- function(image_dirs, dem_path, sitecode, 
                                    output_path=NULL, aoi_file=NULL, n_cpus=1, 
                                    cleartmp=FALSE,  overwrite=FALSE, 
                                    notify=print) {
    if (!file_test("-d", dem_path)) {
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

    for (n in 1:length(image_stacks)) {
        image_stack <- image_stacks[[n]]
        mask_stack <- mask_stacks[[n]]
        band1_imagefile <- image_files[[n]][1]

        if (is.null(output_path)) {
            output_path <- dirname(band1_imagefile)
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
        if (!is.null(aoi_file)) {
            image_basename <- paste(image_basename, 'crop', sep='_')
        }

        timer <- start_timer(timer, label=paste('Preprocessing', image_basename))

        ######################################################################
        # Crop image to landsat path/row, after intersecting it with the 
        # supplied AOI
        timer <- start_timer(timer, label=paste(image_basename, '-', 'crop'))
        if (!is.null(aoi_file)) {
            aoi <- readOGR(dirname(aoi_file), basename(file_path_sans_ext(aoi_file)))
            if (proj4string(aoi) != proj4string(image_stack)) {
                aoi <- spTransform(aoi, CRS(proj4string(image_stack)))
            }
            pathrow_area <- pathrow_poly(as.numeric(WRS_Path), 
                                         as.numeric(WRS_Row))
            if (proj4string(pathrow_area) != proj4string(aoi)) {
                pathrow_area <- spTransform(pathrow_area, 
                                            CRS(proj4string(aoi)))
            }
            crop_area <- gIntersection(pathrow_area, aoi, byid=TRUE)
        } else {
            crop_area <- pathrow_poly(as.numeric(WRS_Path), 
                                      as.numeric(WRS_Row))
            if (proj4string(crop_area) != proj4string(image_stack)) {
                crop_area <- spTransform(crop_area, 
                                         CRS(proj4string(image_stack)))
            }
        }

        image_stack <- crop(image_stack, crop_area)
        image_stack <- mask(image_stack, crop_area)
        image_stack <- extend(image_stack, crop_area)

        mask_stack <- crop(mask_stack, crop_area)
        mask_stack <- mask(mask_stack, crop_area)
        mask_stack <- extend(mask_stack, crop_area)

        timer <- stop_timer(timer, label=paste(image_basename, '-', 'crop'))

        ######################################################################
        # Mask out clouds and missing values
        timer <- start_timer(timer, label=paste(image_basename, '-', 'masking'))

        # The combined cloud mask includes the cloud_QA, cloud_shadow_QA, and 
        # adjacent_cloud_QA layers. Missing or clouded pixels are coded as 0, while 
        # good pixels are coded as 1.
        image_stack_mask_path <- file.path(this_output_path,
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

        image_stack <- mask(image_stack, image_stack_mask, maskvalue=0)

        mask_stack_path <- file.path(this_output_path,
                                     paste(sitecode, image_basename, 
                                           'masks.envi', sep='_'))
        mask_stack <- writeRaster(mask_stack, filename=mask_stack_path, 
                                  overwrite=overwrite, 
                                  datatype=dataType(mask_stack)[1])
        timer <- stop_timer(timer, label=paste(image_basename, '-', 'masking'))

        ######################################################################
        # Load dem, slope, and aspect
        dem_filename <- file.path(dem_path, paste0('dem_', WRS_Path, '-', WRS_Row, 
                                                   '.envi'))
        dem <- raster(dem_filename)

        slopeaspect_filename <- file.path(dem_path,
                                          paste0('slopeaspect_', 
                                                 WRS_Path, '-', WRS_Row, '.envi'))
        slopeaspect <- brick(slopeaspect_filename)

        if (!compareRaster(dem, image_stack, extent=FALSE, rowcol=FALSE, 
                           crs=TRUE, stopiffalse=FALSE)) {
            warning(paste("DEM projection does not match image projection - reprojecting", 
                          image_basename))
            timer <- start_timer(timer, label=paste(image_basename, '-', 'reprojecting DEM'))
            dem <- projectRaster(dem, image_stack)
            timer <- stop_timer(timer, label=paste(image_basename, '-', 'reprojecting DEM'))
            timer <- start_timer(timer, label=paste(image_basename, '-', 'reprojecting slopeaspect'))
            slopeaspect <- projectRaster(slopeaspect, image_stack)
            timer <- stop_timer(timer, label=paste(image_basename, '-', 'reprojecting slopeaspect'))
        }
        # Since the projections match, make sure the proj4strings are identical 
        # so rgeos doesn't throw an error
        proj4string(dem) <- proj4string(image_stack)
        proj4string(slopeaspect) <- proj4string(image_stack)

        ######################################################################
        # Perform topographic correction
        timer <- start_timer(timer, label=paste(image_basename, '-', 'topocorr'))
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
        topocorr_filename <- file.path(this_output_path,
                                       paste(sitecode, image_basename, 
                                             'tc.envi', sep='_'))
        image_stack <- topographic_corr(image_stack, slopeaspect_flt, sunelev, 
                                        sunazimuth, method='minnaert_full', 
                                        asinteger=TRUE, 
                                        sampleindices=sampleindices)
        image_stack <- writeRaster(image_stack, filename=topocorr_filename, 
                                   overwrite=overwrite, datatype='INT2S')
        timer <- stop_timer(timer, label=paste(image_basename, '-', 'topocorr'))

        timer <- stop_timer(timer, label=paste('Preprocessing', image_basename))

        if (cleartmp) removeTmpFiles(h=1)
    }
    if (n_cpus > 1) endCluster()

    timer <- stop_timer(timer, label='Preprocessing images')
}
