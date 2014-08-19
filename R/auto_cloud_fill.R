pct_clouds <- function(cloud_mask) {
    num_clouds <- cellStats(cloud_mask == 1, stat='sum', na.rm=TRUE)
    num_clear <- cellStats(cloud_mask == 0, stat='sum', na.rm=TRUE)
    return((num_clouds / (num_clouds + num_clear)) * 100)
}

#' Automated removal of clouds from Landsat CDR imagery
#'
#' Uses one of four cloud reomval algorithms (see \code{\link{cloud_remove}}) 
#' to remove thick clouds from Landsat imagery. In hilly areas, topographic 
#' correction should be done before cloud fill.
#'
#' The \code{auto_cloud_fill} function allows an analyst to automatically 
#' construct a cloud-filled image after specifying: \code{data_dir} (a folder 
#' of Landsat images), \code{wrspath} and \code{wrsrow} (the WRS-2 path/row to 
#' use), and \code{start_date} and \code{end_date} (a start and end date 
#' limiting the images to use in the algorithm).  The analyst can also 
#' optionally specify a \code{base_date}, and the \code{auto_cloud_fill} 
#' function will automatically pick the image closest to that date to use as 
#' the base image.
#' 
#' As the \code{auto_cloud_fill} function automatically chooses images for 
#' inclusion in the cloud fill process, it relies on having images stored on 
#' disk in a particular way, and currently only supports cloud fill for Landsat 
#' CDR surface reflectance images. To ensure that images are correctly stored 
#' on your hard disk, use the \code{\link{auto_preprocess_landsat}} function to 
#' extract the original Landsat CDR hdf files from the USGS archive. The 
#' \code{auto_preprocess_landsat} function will ensure that images are 
#' extracted and renamed properly so that they can be used with the 
#' \code{auto_cloud_fill} script.
#'
#' @export
#' @importFrom tools file_path_sans_ext
#' @importFrom lubridate as.duration new_interval
#' @importFrom stringr str_extract
#' @importFrom SDMTools ConnCompLabel
#' @param data_dir folder where input images are located, with filenames as 
#' output by the \code{\link{auto_preprocess_landsat}} function. This folder 
#' will be searched recursively for images (taking the below path/row, date, 
#' and topographic correction options into account).
#' @param wrspath World Reference System (WRS) path
#' @param wrsrow World Reference System (WRS) row
#' @param start_date start date of period from which images will be chosen to 
#' fill cloudy areas in the base image (as \code{Date} object)
#' @param end_date end date of period from which images will be chosen to fill 
#' cloudy areas in the the base image (as \code{Date} object)
#' @param base_date ideal date for base image (base image will be chosen as the 
#' image among the available images that is closest to this date). If NULL, 
#' then the base image will be the image with the lowest cloud cover.
#' @param out_name base filename (without an extension - see \code{ext} 
#' argument) for cloud filled image.  The mask file for the cloud filled image 
#' will be saved with the same name, with the added suffix "_mask".
#' @param tc if \code{TRUE}, use topographically corrected imagery as output by 
#' \code{auto_preprocess_landsat}. IF \code{FALSE} use bands 1-5 and 7 surface 
#' reflectance as output by \code{unstack_ledaps} or 
#' \code{auto_preprocess_landsat} (if \code{auto_preprocess_landsat} was also 
#' run with tc=FALSE).
#' @param ext file extension to use when searching for input rasters and when 
#' saving output rasters (determines output file format). Should match file 
#' extension of input rasters (and should most likely match the value chosen 
#' for \code{ext} when \code{auto_preprocess_landsat} was run).
#' @param sensors choose the sensors to include when selecting images (useful 
#' for excluding images from a particular satellite if desired). Can be any of 
#' "L4T", "L5T", "L7E", and/or "L8C".
#' @param img_type type of Landsat imagery to preprocess. Can be "CDR" for 
#' Landsat Climate Data Record (CDR) imagery in HDR format, or "L1T" for 
#' Standard Terrain Correction (Level 1T) imagery. Note that if L1T imagery is 
#' used, fmask must be run locally (see https://code.google.com/p/fmask) prior 
#' to using \code{auto_preprocess_landsat}.
#' @param threshold maximum percent cloud cover allowable in base image. Cloud 
#' fill will iterate until percent cloud cover in base image is below this 
#' value, or until \code{max_iter} iterations have been run
#' @param max_iter maximum number of times to run cloud fill script
#' @param notify notifier to use (defaults to \code{print} function).  See the 
#' \code{notifyR} package for one way of sending notifications from R.  The 
#' \code{notify} function should accept a string as the only argument.
#' @param verbose whether to print detailed status messages. Set to FALSE or 0 
#' for no status messages. Set to 1 for basic status messages. Set to 2 for 
#' detailed status messages.
#' @param overwrite whether to overwrite \code{out_name} if it already exists
#' @param ...  additional arguments passed to \code{\link{cloud_remove}}, such 
#' as \code{DN_min}, \code{DN_max}, \code{algorithm}, \code{byblock}, 
#' \code{verbose}, etc. See \code{\link{cloud_remove}} for details
#' @return a list with two elements: "filled", a \code{Raster*} object with 
#' cloud filled image, and "mask", a \code{RasterLayer} object with the cloud 
#' mask for the cloud filled image.
#' @references Zhu, X., Gao, F., Liu, D., Chen, J., 2012. A modified 
#' neighborhood similar pixel interpolator approach for removing thick clouds 
#' in Landsat images.  Geoscience and Remote Sensing Letters, IEEE 9, 521--525.  
#' doi:10.1109/LGRS.2011.2173290
auto_cloud_fill <- function(data_dir, wrspath, wrsrow, start_date, end_date, 
                            out_name, base_date=NULL, tc=TRUE, ext='tif',
                            sensors=c('L4T', 'L5T', 'L7E', 'L8C'), 
                            img_type="CDR", threshold=1, max_iter=5, 
                            notify=print, verbose=1, overwrite=FALSE, ...) {
    if (!file_test('-d', data_dir)) {
        stop('data_dir does not exist')
    }
    if (!file_test('-d', dirname(out_name))) {
        stop('output folder does not exist')
    }
    if (file_path_sans_ext(out_name) != out_name) {
        stop('out_name should not have a file extension')
    }

    ext <- gsub('^[.]', '', ext)

    output_file <- paste0(out_name, '.', ext)
    if (file_test('-f', output_file) & !overwrite) {
        stop(paste0('output file "', output_file, '" already exists'))
    }
    if (!all(sensors %in% c('L4T', 'L5T', 'L7E', 'L8C'))) {
        stop('"sensors" must be a list of one or more of: "L4T", "L5T", "L7E", "L8C"')
    }

    log_file <- file(paste0(out_name, '_log.txt'), open="wt")
    msg <- function(txt) {
        cat(paste0(txt, '\n'), file=log_file, append=TRUE)
        print(txt)
    }

    timer <- Track_time(msg)
    timer <- start_timer(timer, label='Cloud fill')

    stopifnot(class(start_date) == 'Date')
    stopifnot(class(end_date) == 'Date')

    wrspath <- sprintf('%03i', wrspath)
    wrsrow <- sprintf('%03i', wrsrow)

    # Find image files based on start and end dates
    prefix_re <- "^([a-zA-Z]*_)?"
    #pathrow_re <-"[012][0-9]{2}-[012][0-9]{2}"
    pathrow_re <- paste(wrspath, wrsrow, sep='-')
    date_re <-"((19)|(2[01]))[0-9]{2}-[0123][0-9]{2}"
    if (img_type == "CDR") {
        sensor_re <- paste0('(', paste0(paste0('(', sensors,')'), collapse='|'), ')', "SR")
    } else if (img_type == "L1T") {
        sensor_re <- paste0('(', paste0(paste0('(', sensors,')'), collapse='|'), ')', "L1T")
    } else {
        stop(paste(img_type, "is not a recognized img_type"))
    }
    if (tc) {
        suffix_re <- paste0('_tc.', ext, '$')
    } else {
        suffix_re <- paste0('.', ext, '$')
    }
    file_re <- paste0(prefix_re, paste(pathrow_re, date_re, sensor_re, 
                                       sep='_'), suffix_re)
    img_files <- dir(data_dir, pattern=file_re, recursive=TRUE)

    img_dates <- str_extract(basename(img_files), date_re)
    img_dates <- as.Date(img_dates, '%Y-%j')

    which_files <- which((img_dates >= start_date) &
                          (img_dates < end_date))
    img_dates <- img_dates[which_files]
    img_files <- file.path(data_dir, img_files[which_files])

    if (length(img_files) == 0) {
        stop('no images found - check date_dir, check wrspath, wrsrow, start_date, and end_date')
    } else if (length(img_files) < 2) {
        stop(paste('Only', length(img_files),
                   'image(s) found. Need at least two images to perform cloud fill'))
    }

    if (verbose > 0) {
        msg(paste('Found', length(img_files), 'image(s)'))
        timer <- start_timer(timer, label='Analyzing cloud cover in input images')
    }
    # Run QA stats
    fmasks <- list()
    fill_QAs <- list()
    imgs <- list()
    for (img_file in img_files) {
        masks_file <- paste0(file_path_sans_ext(img_file), '_masks.', ext)
        if (!file_test('-f', masks_file)) {
            masks_file <- gsub(suffix_re, paste0('_masks.', ext), img_file)
            if (file_test('-f', masks_file)) {
                warning('using masks file with old format (pre v0.5) teamlucc naming')
            } else {
                stop('could not find masks file')
            }
        }
        this_fill_QA <- raster(masks_file, band=1)
        fill_QAs <- c(fill_QAs, this_fill_QA)
        this_fmask <- raster(masks_file, band=2)
        fmasks <- c(fmasks, this_fmask)
        this_img <- stack(img_file)
        imgs <- c(imgs, stack(this_img))
    }
    
    compareRaster(imgs, res=TRUE, orig=TRUE)
    compareRaster(fmasks, res=TRUE, orig=TRUE)

    freq_table <- freq(stack(fmasks), useNA='no', merge=TRUE)
    # Convert frequency table to fractions
    freq_table[-1] <- freq_table[-1] / colSums(freq_table[-1], na.rm=TRUE)
    if (verbose > 0) {
        timer <- stop_timer(timer, label='Analyzing cloud cover in input images')
    }
    if (verbose > 0) {
        timer <- start_timer(timer, label='Calculating cloud masks')
    }

    # Find image that is either closest to base date, or has the maximum 
    # percent clear
    if (is.null(base_date)) {
        clear_row <- which(freq_table$value == 0)
        base_img_index <- which(freq_table[clear_row, -1] == 
                                max(freq_table[clear_row, -1]))
    } else {
        base_date_diff <- lapply(img_dates, function(x) 
                                 as.duration(new_interval(x, base_date)))
        base_date_diff <- abs(unlist(base_date_diff))
        base_img_index <- which(base_date_diff == min(base_date_diff))
        # Handle ties - two images that are the same distance from base date.  
        # Default to earlier image.
        if (length(base_img_index) > 1) {
            base_img_index <- base_img_index[1]
        }
    }

    # Save the original base image fmask so it can be used to recode the final 
    # cloud mask at the end of cloud filling
    base_fmask <- fmasks[[base_img_index]]
    base_fill_QA <- fill_QAs[[base_img_index]]

    # Convert masks to indicate: 0 = clear; 1 = cloud or shadow; 2 = fill
    #
    #   fmask_band key:
    #       0 = clear
    #       1 = water
    #       2 = cloud_shadow
    #       3 = snow
    #       4 = cloud
    #       255 = fill value
    calc_cloud_mask <- function(fmask, img) {
        # Code clouds and cloud shadows as 1
        ret <- (fmask == 2) | (fmask == 4)
        # Code fill as 2
        ret[fmask == 255] <- 2
        # Code other missing data that is not in fill areas as NA. The (ret != 
        # 1) test is necessary to ensures that only NAs that are NOT in clouds 
        # will be copied to the mask images (the assumption being that NAs in 
        # clouds should be marked as cloud and fill should be attempted).  This 
        # is necessary in case clouded areas in img are mistakenly coded NA 
        # (they should not be).
        ret[(ret != 1) & (ret != 2) & is.na(img)] <- NA
        return(ret)
    }
    for (n in 1:length(fmasks)) {
        fmasks[n] <- overlay(fmasks[[n]], imgs[[n]][[1]], fun=calc_cloud_mask, 
                             datatype=dataType(fmasks[[n]]))
    }

    base_img <- imgs[[base_img_index]]
    imgs <- imgs[-base_img_index]
    base_mask <- fmasks[[base_img_index]]
    fmasks <- fmasks[-base_img_index]

    base_img_date <- img_dates[base_img_index]
    img_dates <- img_dates[-base_img_index]

    if (verbose > 0) {
        msg(paste('Using image from', base_img_date, 'as base image.'))
    }

    if (verbose > 0) {
        timer <- stop_timer(timer, label='Calculating cloud masks')
    }

    if (verbose > 0) {
        timer <- start_timer(timer, label='Masking base image')
    }
    # Mask out clouds in base image. Save this image to disk so it is available 
    # even if no cloud fill is done (if the pct_clouds in this image is below 
    # the threshold).
    base_img <- overlay(base_img, base_mask,
        fun=function(base_vals, mask_vals) {
            # Set clouds/shadows to 0
            base_vals[mask_vals == 1] <- 0
            # Allow fill to be attempted in NA areas
            base_vals[is.na(base_vals)] <- 0
            # Set slc-off gaps and areas outside scene to NA
            base_vals[mask_vals == 2] <- NA
            return(base_vals)
        }, datatype=dataType(base_img[[1]]), 
        filename=extension(rasterTmpFile(), ext), overwrite=overwrite)

    cur_pct_clouds <- pct_clouds(base_mask)

    if (verbose > 0) {
        msg(paste0('Base image has ', round(cur_pct_clouds, 2), '% cloud cover before fill'))
    }

    if (verbose > 0) {
        timer <- stop_timer(timer, label='Masking base image')
    }

    n <- 0
    while ((cur_pct_clouds > threshold) & (n < max_iter) & (length(imgs) >= 1)) {
        if (verbose > 0) {
            timer <- start_timer(timer, label=paste('Fill iteration', n + 1))
        }

        # Calculate a raster indicating the pixels in each potential fill image 
        # that are available for filling pixels of base_img that are missing 
        # due to cloud contamination. Areas coded 1 are missing due to cloud or 
        # shadow in the base image and are available in the merge image. This 
        # will return a stack with number of layers equal to number of masks.
        fill_areas <- overlay(base_mask, stack(fmasks),
            fun=function(base_mask_vals, fill_mask_vals) {
                ret <- rep(NA, length(base_mask_vals))
                # Code cloudy in base, clear in fill as 1
                ret[(base_mask_vals == 1) & (fill_mask_vals == 0)] <- 1
                # Code clear in base, clear in fill as 0
                ret[(base_mask_vals == 0) & (fill_mask_vals == 0)] <- 0
                # Code NA in base, clear in fill as clouded, so these NAs will 
                # be filled if possible.
                ret[is.na(base_mask_vals) & (fill_mask_vals == 0)] <- 1
                # Ensure SLC-off gaps and background areas in each image are 
                # not filled:
                ret[(base_mask_vals == 2) | (fill_mask_vals == 2)] <- NA
                return(ret)
            }, datatype=dataType(base_mask))
        fill_areas_freq <- freq(fill_areas, useNA='no', merge=TRUE)
        # Below is necessary as for some reason when fill_areas is of length 
        # one, freq returns a matrix rather than a data.frame
        fill_areas_freq <- as.data.frame(fill_areas_freq)
        # Select the fill image with the maximum number of available pixels 
        # (counting only pixels in the fill image that are not ALSO clouded in 
        # the fill image)
        avail_fill_row <- which(fill_areas_freq$value == 1)
        if (length(avail_fill_row) == 0) {
            msg(paste('No fill pixels available. Stopping fill.'))
            break
        }
        # Remove the now unnecessary "value" column
        fill_areas_freq <- fill_areas_freq[!(names(fill_areas_freq) == 'value')]
        fill_img_index <- which(fill_areas_freq[avail_fill_row, ] == 
                                max(fill_areas_freq[avail_fill_row, ], na.rm=TRUE))
        if ((length(fill_img_index) == 0) ||
            (fill_areas_freq[avail_fill_row, fill_img_index] == 0)) {
            msg(paste('No fill pixels available. Stopping fill.'))
            break
        }

        fill_img <- imgs[[fill_img_index]]
        imgs <- imgs[-fill_img_index]
        base_img_mask <- fill_areas[[fill_img_index]]
        fmasks <- fmasks[-fill_img_index]
        fill_img_date <- img_dates[fill_img_index]
        img_dates <- img_dates[-fill_img_index]

        # Add numbered IDs to the cloud patches
        base_img_mask <- ConnCompLabel(base_img_mask)

        # Ensure dataType is properly set prior to handing off to IDL
        dataType(base_img_mask) <- 'INT2S'

        if (verbose > 0) {
            msg(paste0('Filling image from ', base_img_date,
                          ' with image from ', fill_img_date, '.'))
            timer <- start_timer(timer, label="Performing fill")
        }
        base_img <- cloud_remove(base_img, fill_img, base_img_mask, 
                                 out_name=extension(rasterTmpFile(), ext), 
                                 verbose=verbose, overwrite=TRUE, ...)
        # base_img <- cloud_remove(base_img, fill_img, base_img_mask, 
        #                          out_name=extension(rasterTmpFile(), ext), 
        #                          verbose=verbose, overwrite=TRUE, 
        #                          DN_min=DN_min, DN_max=DN_max, 
        #                          algorithm=algorithm, byblock=byblock)
        if (verbose > 0) {
            timer <- stop_timer(timer, label="Performing fill")
        }

        # Revise base mask to account for newly filled pixels
        base_mask <- overlay(base_mask, base_img[[1]],
            fun=function(mask_vals, filled_vals) {
                mask_vals[(mask_vals == 1) & (filled_vals != 0)] <- 0
                return(mask_vals)
            }, datatype=dataType(base_mask), 
            filename=extension(rasterTmpFile(), ext), overwrite=TRUE)

        cur_pct_clouds <- pct_clouds(base_mask)
        if (verbose > 0) {
            msg(paste0('Base image has ', round(cur_pct_clouds, 2),
                          '% cloud cover remaining'))
            timer <- stop_timer(timer, label=paste('Fill iteration', n + 1))
        }

        n <- n + 1
    }

    base_img <- writeRaster(base_img, filename=output_file, datatype="INT2S", 
                            overwrite=overwrite)

    # Recode base mask so final coding matches that of fmask (though cloud and 
    # cloud shadow are no longer differentiated)
    #   fmask_band key:
    #       0 = clear
    #       1 = water
    #       2 = cloud_shadow
    #       3 = snow
    #       4 = cloud
    #       255 = fill value
    #   base_mask key:
    #   	0 = clear
    #   	1 = cloud
    #   	2 = fill
    mask_output_file <- paste0(out_name, '_masks.', ext)
    filled_fmask <- overlay(base_mask, base_fmask,
        fun=function(after_fill, before_fill) {
            ret <- after_fill
            # Code clear after filling but water in fmask as water (1 in fmask)
            ret[(after_fill == 0) & (before_fill == 1)] <- 1
            # Code clear after filling but snow in fmask as snow (3 in fmask)
            ret[(after_fill == 0) & (before_fill == 3)] <- 3
            # Code cloudy after filling as cloud (4 in fmask)
            ret[after_fill == 1] <- 4
            # Code gap in fmask as gap (3 in fmask)
            ret[before_fill == 255] <- 255
            return(ret)
        }, datatype=dataType(base_mask))
    final_masks <- stack(base_fill_QA, filled_fmask)
    names(final_masks) <- c("fill_QA", "fmask")
    final_masks <- writeRaster(final_masks, datatype=dataType(base_mask), 
                               filename=mask_output_file, overwrite=TRUE)

    timer <- stop_timer(timer, label='Cloud fill')

    close(log_file)

    return(list(filled=base_img, mask=final_masks))
}
