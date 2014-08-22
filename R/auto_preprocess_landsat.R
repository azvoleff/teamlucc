get_gdalinfo_item <- function(item, gdalinfo_text) {
    gdalinfo_text <- gdalinfo_text[grepl(paste0('^[ ]*', item), gdalinfo_text)]
    if (length(gdalinfo_text) > 1) stop('more than one item found')
    gdalinfo_text <- gsub(paste0('[ ]*', item, '='), '', gdalinfo_text)
    return(gdalinfo_text)
}

get_mtl_item <- function(item, mtl_txt) {
    mtl_txt <- mtl_txt[grepl(paste0('^[ ]*', item), mtl_txt)]
    if (length(mtl_txt) > 1) stop('more than one item found')
    mtl_txt <- gsub(paste0('[ ]*', item, ' = '), '', mtl_txt)
    # Strip leading/following quotes
    mtl_txt <- gsub('^"', '', mtl_txt)
    mtl_txt <- gsub('"$', '', mtl_txt)
    return(mtl_txt)
}

#' @importFrom stringr str_extract
#' @importFrom gdalUtils gdalinfo
get_metadata <- function(ls_file, img_type) {
    meta <- list()
    if (img_type == "CDR") {
        ls_file_gdalinfo <- gdalinfo(ls_file)
        aq_date <- get_gdalinfo_item('AcquisitionDate', ls_file_gdalinfo)
        meta$aq_date <- strptime(aq_date, format="%Y-%m-%dT%H:%M:%OSZ", tz="UTC")
        meta$WRS_Path <- sprintf('%03i', as.numeric(get_gdalinfo_item('WRS_Path', ls_file_gdalinfo)))
        meta$WRS_Row <- sprintf('%03i', as.numeric(get_gdalinfo_item('WRS_Row', ls_file_gdalinfo)))
        meta$sunelev <- 90 - as.numeric(get_gdalinfo_item('SolarZenith', ls_file_gdalinfo))
        meta$sunazimuth <- as.numeric(get_gdalinfo_item('SolarAzimuth', ls_file_gdalinfo))
        meta$short_name  <- get_gdalinfo_item('ShortName', ls_file_gdalinfo)
    } else if (img_type == "L1T") {
        if (!grepl("_MTL.txt$", ls_file)) {
            stop("ls_file must be a *_MTL.txt file")
        }
        mtl_txt <- readLines(ls_file, warn=FALSE)
        aq_date <- get_mtl_item('DATE_ACQUIRED', mtl_txt)
        aq_time <- get_mtl_item('SCENE_CENTER_TIME', mtl_txt)
        meta$aq_date <- strptime(paste0(aq_date, "T", aq_time), format="%Y-%m-%dT%H:%M:%OSZ", tz="UTC")
        meta$WRS_Path <- sprintf('%03i', as.numeric(get_mtl_item('WRS_PATH', mtl_txt)))
        meta$WRS_Row <- sprintf('%03i', as.numeric(get_mtl_item('WRS_ROW', mtl_txt)))
        meta$sunelev <- as.numeric(get_mtl_item('SUN_ELEVATION', mtl_txt))
        meta$sunazimuth <- as.numeric(get_mtl_item('SUN_AZIMUTH', mtl_txt))
        # Build a shortname based on satellite and img_type that is consistent 
        # with the format of the CDR image shortnames
        satellite <- str_extract(get_mtl_item('SPACECRAFT_ID', mtl_txt), '[4578]')
        sensor_string <- str_extract(basename(ls_file), '^((LT[45])|(LE7)|(LC8))')
        meta$short_name  <- paste0(substr(sensor_string, 1, 1),
                                   substr(sensor_string, 3, 3),
                                   substr(sensor_string, 2, 2), img_type)
    } else {
        stop(paste(img_type, "is not a recognized img_type"))
    }
    return(meta)
}

calc_cloud_mask <- function(mask_stack, mask_type, ...) {
    if (mask_type == 'fmask') {
        # Make a mask where clouds and gaps are coded as 1, clear as 0
        # fmask_band key:
        # 	0 = clear
        # 	1 = water
        # 	2 = cloud_shadow
        # 	3 = snow
        # 	4 = cloud
        # 	255 = fill value
        cloud_mask <- calc(mask_stack$fmask_band,
            fun=function(fmask) {
                return((fmask == 2) | (fmask == 4) | (fmask == 255))
            }, datatype='INT2S', ...)
    } else if (mask_type == '6S') {
        # This cloud mask includes the cloud_QA, cloud_shadow_QA, and 
        # adjacent_cloud_QA layers. Pixels in cloud, cloud shadow, or 
        # adjacent cloud are coded as 1.
        cloud_mask <- overlay(mask_stack$fill_QA,
                              mask_stack$cloud_QA, 
                              mask_stack$cloud_shadow_QA, 
                              mask_stack$adjacent_cloud_QA,
            fun=function(fill, clo, sha, adj) {
                return((fill == 255) | (clo == 255) | (sha == 255) | 
                       (adj == 255))
            }, datatype='INT2S', ...)

    } else if (mask_type == 'both') {
        cloud_mask <- overlay(mask_stack$fmask_band, 
                              mask_stack$cloud_QA, 
                              mask_stack$cloud_shadow_QA, 
                              mask_stack$adjacent_cloud_QA,
            fun=function(fmask, clo, sha, adj) {
                return((fmask == 2) | (fmask == 4) | (fmask == 255) | 
                       (clo == 255) | (sha == 255) | (adj == 255))
            }, datatype='INT2S', ...)
    } else {
        stop(paste0('unrecognized option "', cloud_mask, '" for mask_type"'))
    }
    return(cloud_mask)
}

#' @importFrom gdalUtils get_subdatasets gdalbuildvrt
build_band_vrt <- function(ls_file, band_vrt_file, img_type) {
    image_bands <- c('band1', 'band2', 'band3', 'band4', 'band5', 'band7')
    if (img_type == "CDR") {
        sds <- get_subdatasets(ls_file)
        band_sds <- sds[grepl(paste0(':(', paste(image_bands, collapse='|'), ')$'), sds)]
        gdalbuildvrt(band_sds, band_vrt_file, separate=TRUE)
    } else if (img_type == "L1T") {
        if (!grepl("_MTL.txt$", ls_file)) {
            stop("ls_file must be a *_MTL.txt file")
        }
        ls_file_base <- gsub("_MTL.txt", "", ls_file)
        ls_files <- dir(dirname(ls_file_base),
                        pattern=paste0(basename(ls_file_base), '_B[123457].((TIF)|(tif))$'),
                        full.names=TRUE)
        gdalbuildvrt(ls_files, band_vrt_file, separate=TRUE)

    } else {
        stop(paste(img_type, "is not a recognized img_type"))
    }
    return(image_bands)
}

#' @importFrom gdalUtils get_subdatasets gdalbuildvrt
build_mask_vrt <- function(ls_file, mask_vrt_file, img_type) {
    if (img_type == "CDR") {
        mask_bands <- c('fill_QA', 'cfmask_band', 'cloud_QA', 'cloud_shadow_QA', 
                        'adjacent_cloud_QA')
        sds <- get_subdatasets(ls_file)
        # Below is to support CDR imagery downloaded prior to late August 2014
        if (any(grepl("fmask_band", sds))) {
            warning('Using "fmask_band" instead of newer "cfmask_band" band name')
            mask_bands[grepl("^cfmask_band$", mask_bands)] <- "fmask_band"
        }
        mask_sds <- sds[grepl(paste0(':(', paste(mask_bands, collapse='|'), ')$'), sds)]
        stopifnot(length(mask_sds) == 5)
        gdalbuildvrt(mask_sds, mask_vrt_file, separate=TRUE, srcnodata='None')
    } else if (img_type == "L1T") {
        mask_bands <- c('fill_QA', 'fmask_band')
        if (!grepl("_MTL.txt$", ls_file)) {
            stop("ls_file must be a *_MTL.txt file")
        }
        ls_file_base <- gsub("_MTL.txt", "", ls_file)

        fmask_file <- dir(dirname(ls_file_base),
                          pattern=paste0(basename(ls_file_base), '_MTLFmask$'),
                          full.names=TRUE)

        # Calculate a QA mask file from the fmask file, since teamlucc expects 
        # this file as part of the mask stack.
        qa_mask_file <- extension(rasterTmpFile(), '.tif')
        # TODO: Check if this is proper coding - should it be reversed?
        qa_mask <- calc(raster(fmask_file),
                        fun=function(x) {
                            out <- x == 255
                            out[x == 255] <- 255
                            return(out)
                        }, datatype="INT2S", filename=qa_mask_file)

        # Note that allow_projection_difference is used below as GDAL thinks 
        # the two images have different projection systems, even though they 
        # are in identical projection systems.
        gdalbuildvrt(c(qa_mask_file, fmask_file), mask_vrt_file, 
                     separate=TRUE, allow_projection_difference=TRUE,
                     srcnodata='None')
    } else {
        stop(paste(img_type, "is not a recognized img_type"))
    }
    return(mask_bands)
}

#' Preprocess surface reflectance imagery from the Landsat CDR archive
#'
#' This function preprocesses surface reflectance imagery from the Landsat 
#' Climate Data Record (CDR) archive. \code{auto_preprocess_landsat} can 
#' reproject CDR tiles to match the projection of a given \code{aoi}, crop the 
#' tiles to match the \code{aoi} or a common WRS-2 path/row polygon, mask 
#' missing data and clouds out of the CDR tiles, and perform topographic 
#' correction.
#'
#' \code{mask_type} chooses the cloud mask to use if topographic correction is 
#' performed (\code{tc=TRUE}). The mask can be one of three different options: 
#' "6S", "fmask", or "combined". Each option uses a different combination of 
#' cloud mask layers from the CDR product. The "6S" masks out any areas coded 
#' as fill (fill_QA=255), cloud (cloud_QA=255), cloud shadow
#' (cloud_shadow_QA=255) or adjacent to cloud (adjacent_cloud_QA=255). The 
#' "fmask" option masks out any areas coded as fill (fmask=255), cloud 
#' (fmask=4) or cloud shadow (fmask=2).  The combined option combines the "6S" 
#' and "fmask" approaches to masks out areas coded as fill, cloud, cloud 
#' shadow, or adjacent to cloud using either method. Note that "fmask" is the 
#' only supported option when \code{img_type} is L1T.
#'
#' Prior to running \code{auto_preprocess_landsat}, \code{\link{espa_extract}} 
#' should be used to extract the original zipfiles supplied by USGS. To perform 
#' topographic correction with \code{auto_preprocess_landsat}, first run 
#' \code{\link{auto_setup_dem}} to preprocess a set of DEM tiles. Then run 
#' \code{auto_preprocess_landsat} with the \code{tc=TRUE} option.
#'
#' If topographic correction is being performed, it will be run in parallel if 
#' a parallel backend is registered with \code{\link{foreach}}.
#'
#' @export
#' @importFrom rgeos gIntersection
#' @importFrom wrspathrow pathrow_poly
#' @importFrom tools file_path_sans_ext
#' @importFrom gdalUtils gdalwarp
#' @importFrom sp is.projected
#' @param image_dirs list of paths to a set of Landsat CDR image files in HDF 
#' format
#' @param prefix string to use as a prefix for all filenames
#' @param img_type type of Landsat imagery to preprocess. Can be "CDR" for 
#' Landsat Climate Data Record (CDR) imagery in HDR format, or "L1T" for 
#' Standard Terrain Correction (Level 1T) imagery. Note that if L1T imagery is 
#' used, fmask must be run locally (see https://code.google.com/p/fmask) prior 
#' to using \code{auto_preprocess_landsat}.
#' @param tc whether to topographically correct imagery (if \code{TRUE}, then 
#' \code{dem_path} must be specified)
#' @param dem_path path to a set of DEMs as output by \code{auto_setup_dem} 
#' (only required if tc=TRUE)
#' @param aoi area of interest (AOI), as a \code{SpatialPolygonsDataFrame}.  If 
#' supplied, this aoi is used to crop and set the projection system of the 
#' output. Must be in a projected coordinate system.
#' @param output_path the path to use for the output (optional - if NULL then 
#' output images will be saved alongside the input images in the same folder).
#' @param mask_type which cloud mask to use to mask clouds when performing 
#' topographic correction. Can be one of "fmask", "6S", or "both".  See 
#' Details.  (Ignored if \code{tc=FALSE)}.
#' @param mask_output if \code{TRUE}, cloud, cloud shadow, and fill areas 
#' (SLC-off gaps and areas with no data) will be set to \code{NA} in the 
#' output. Note this setting affects the final output file only - cloud, cloud 
#' shadow, and gap areas are masked out of the image during topographic 
#' correction regardless of the value of \code{mask_output}.
#' @param n_cpus the number of CPUs to use for processes that can run in 
#' parallel
#' @param cleartmp whether to clear temp files on each run through the loop
#' @param overwrite whether to overwrite existing files (otherwise an error 
#' will be raised)
#' @param of output format to use when saving output rasters. See description 
#' of \code{of} in \code{\link{gdalwarp}}.
#' @param ext file extension to use when saving output rasters (determines 
#' output file format). Should match file extension for output format chosen by 
#' \code{of}.
#' @param notify notifier to use (defaults to \code{print} function). See the 
#' \code{notifyR} package for one way of sending notifications from R. The 
#' \code{notify} function should accept a string as the only argument.
#' @param verbose whether to print detailed status messages and timing 
#' information
#' @return nothing - used for the side effect of preprocessing imagery
#' @seealso \code{\link{espa_extract}}, \code{\link{unstack_ledapscdr}}, 
#' \code{\link{auto_setup_dem}}
auto_preprocess_landsat <- function(image_dirs, prefix, img_type="CDR", 
                                    tc=FALSE, dem_path=NULL, aoi=NULL, 
                                    output_path=NULL, mask_type='fmask', 
                                    mask_output=FALSE, n_cpus=1, 
                                    cleartmp=FALSE,  overwrite=FALSE, 
                                    of="GTiff", ext='tif', notify=print, 
                                    verbose=FALSE) {
    if (grepl('_', prefix)) {
        stop('prefix cannot contain underscores (_)')
    }
    if (tc && is.null(dem_path)) {
        stop('dem_path must be supplied if tc=TRUE')
    }
    if (tc && !file_test("-d", dem_path)) {
        stop(paste(dem_path, "does not exist"))
    }
    if (!is.null(output_path) && !file_test("-d", output_path)) {
        stop(paste(output_path, "does not exist"))
    }
    if (!is.null(aoi)) {
        if (length(aoi) > 1) {
            stop('aoi should be a SpatialPolygonsDataFrame of length 1')
        }
        stopifnot(is.projected(aoi))
    }

    ext <- gsub('^[.]', '', ext)

    # Setup a regex to identify Landsat CDR images
    if (img_type == "CDR") {
        ls_regex <- '^(lndsr.)?((LT4)|(LT5)|(LE7)|(LC8))[0-9]{6}[12][0-9]{6}[a-zA-Z]{3}[0-9]{2}.hdf$'
    } else if (img_type == "L1T") {
        ls_regex <- '((LT[45])|(LE7)|(LC8))[0-9]{6}[12][0-9]{6}[a-zA-Z]{3}[0-9]{2}_MTL.txt$'
    } else {
        stop(paste(img_type, "is not a recognized img_type"))
    }

    if (img_type == "CDR") {
        stopifnot(mask_type %in% c('fmask', '6S', 'both'))
    } else if (img_type == "L1T") {
        stopifnot(mask_type == 'fmask')
    }

    ls_files <- c()
    for (image_dir in image_dirs) {
        if (!file_test("-d", image_dir)) {
            stop(paste(image_dir, "does not exist"))
        }
        ls_files <- c(ls_files, dir(image_dir, pattern=ls_regex, full.names=TRUE))
    }

    if (length(ls_files) == 0) {
        stop(paste0('No Landsat files found using img_type="', img_type, '".'))
    }

    for (ls_file in ls_files) {
        ######################################################################
        # Determine image basename for use in naming subsequent files
        meta <- get_metadata(ls_file, img_type)

        image_basename <- paste0(meta$WRS_Path, '-', meta$WRS_Row, '_',
                                 format(meta$aq_date, '%Y-%j'), '_', meta$short_name)

        if (is.null(output_path)) {
            this_output_path <- dirname(ls_file)
        } else {
            this_output_path  <- output_path
        }

        if (tc) {
            output_filename <- file.path(this_output_path,
                                         paste0(prefix, '_', image_basename, 
                                                '_tc.', ext))
        } else {
            # Skip topographic correction, so don't append _tc to filename
            output_filename <- file.path(this_output_path,
                                         paste0(prefix, '_', image_basename, 
                                                '.', ext))
        }

        log_file <- file(paste0(file_path_sans_ext(output_filename), '_log.txt'), open="wt")
        msg <- function(txt) {
            cat(paste0(txt, '\n'), file=log_file, append=TRUE)
            print(txt)
        }

        timer <- Track_time(msg)
        timer <- start_timer(timer, label=paste('Preprocessing', image_basename))

        #######################################################################
        # Crop and reproject images to match the projection being used for this 
        # image.  This is either the projection of the aoi (if aoi is 
        # supplied), or the UTM zone of the centroid of this path and row.
        if (verbose) timer <- start_timer(timer, label='cropping and reprojecting')

        band_vrt_file <- extension(rasterTmpFile(), '.vrt')
        band_names <- build_band_vrt(ls_file, band_vrt_file, img_type)
        mask_vrt_file <- extension(rasterTmpFile(), '.vrt')
        mask_band_names <- build_mask_vrt(ls_file, mask_vrt_file, img_type)

        this_pathrow_poly <- pathrow_poly(as.numeric(meta$WRS_Path), 
                                          as.numeric(meta$WRS_Row))
        if (!is.null(aoi)) {
            to_srs <- proj4string(aoi)
        } else {
            to_srs <- utm_zone(this_pathrow_poly, proj4string=TRUE)
        }

        # Calculate minimum bounding box coordinates:
        this_pathrow_poly <- spTransform(this_pathrow_poly, CRS(to_srs))
        if (!is.null(aoi)) {
            # If an aoi IS supplied, match the image extent to that of the AOI 
            # cropped to the appropriate Landsat path/row polygon.
            crop_area <- gIntersection(this_pathrow_poly, aoi, byid=TRUE)
        } else {
            # If an aoi IS NOT supplied, match the image extent to the 
            # appropriate Landsat path/row polygon.
            crop_area <- this_pathrow_poly
        }
        out_te <- as.numeric(bbox(crop_area))

        # Ensure origin is set at 0,0
        to_res <- c(30, 30)
        out_te <- normalize_extent(out_te, to_res)

        image_stack_reproj_file <- extension(rasterTmpFile(), ext)
        image_stack <- gdalwarp(band_vrt_file,
                                dstfile=image_stack_reproj_file,
                                te=out_te, t_srs=to_srs, tr=to_res, 
                                r='cubicspline', output_Raster=TRUE, of=of, 
                                multi=TRUE, wo=paste0("NUM_THREADS=", n_cpus), 
                                overwrite=overwrite, ot='Int16')
        names(image_stack) <- band_names

        mask_stack_reproj_file <- extension(rasterTmpFile(), paste0('.', ext))
        mask_stack <- gdalwarp(mask_vrt_file,
                               dstfile=mask_stack_reproj_file,
                               te=out_te, t_srs=to_srs, tr=to_res, 
                               r='near', output_Raster=TRUE, of=of, 
                               multi=TRUE, wo=paste0("NUM_THREADS=", n_cpus), 
                               overwrite=overwrite, ot='Int16')
        # Can't just directly assign mask_bands as the names since the bands 
        # may have been read in different order from the HDF file
        names(mask_stack) <- mask_band_names

        if (verbose) timer <- stop_timer(timer, label='cropping and reprojecting')

        ######################################################################
        # Perform topographic correction if tc=TRUE
        if (tc) {
            if (verbose) timer <- start_timer(timer, label='topocorr')

            ######################################################################
            # Load dem, slope, and aspect
            slopeaspect_filename <- file.path(dem_path,
                                              paste0('slopeaspect_', 
                                                     meta$WRS_Path, '-', meta$WRS_Row, '.', ext))
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

            compareRaster(slopeaspect, image_stack, orig=TRUE)

            image_stack_mask <- calc_cloud_mask(mask_stack, mask_type)

            image_stack_masked <- image_stack
            image_stack_masked[image_stack_mask] <- NA
            if (ncell(image_stack_masked) > 500000) {
                # Draw a sample for the Minnaert k regression. Note that 
                # sampleRegular with cells=TRUE returns cell numbers in the 
                # first column
                sampleindices <- sampleRegular(image_stack_masked, size=500000, 
                                               cells=TRUE)
                sampleindices <- as.vector(sampleindices[, 1])
            } else {
                sampleindices <- NULL
            }
            # Remember that slopeaspect layers are scaled to INT2S, but 
            # topographic_corr expects them as floats, so apply the scale factors 
            # used in auto_setup_dem
            slopeaspect_flt <- stack(raster(slopeaspect, layer=1) / 10000,
                                     raster(slopeaspect, layer=2) / 1000)
            image_stack_tc <- topographic_corr(image_stack_masked, 
                                               slopeaspect_flt, meta$sunelev, 
                                               meta$sunazimuth, 
                                               method='minnaert_full', 
                                               asinteger=TRUE, 
                                               sampleindices=sampleindices)
            if (!mask_output) {
                # Add back in the original values of areas that were masked out 
                # from the topographic correction:
                image_stack_tc[image_stack_mask] <- image_stack[image_stack_mask]
            }
            image_stack <- image_stack_tc
            
            if (verbose) timer <- stop_timer(timer, label='topocorr')
        }

        ######################################################################
        # Write final data
        if (verbose) timer <- start_timer(timer, label='writing data')

        mask_stack_path <- paste0(file_path_sans_ext(output_filename), 
                                  '_masks.', ext)
        mask_stack <- writeRaster(stack(mask_stack$fill_QA,
                                        mask_stack$fmask_band),
                                  filename=mask_stack_path, 
                                  overwrite=overwrite, datatype='INT2S')
        names(mask_stack) <- c('fill_QA', 'fmask_band')

        image_stack <- writeRaster(image_stack, filename=output_filename, 
                                   overwrite=overwrite, datatype='INT2S')
        if (verbose) timer <- stop_timer(timer, label='writing data')

        timer <- stop_timer(timer, label=paste('Preprocessing', image_basename))

        close(log_file)

        if (cleartmp) removeTmpFiles(h=1)
    }
}
