#' Perform cloud fill for topographically corrected Landsat imagery
#'
#' Uses the GNSPI algorithm from Xioalin Zhu. See \code{\link{fill_clouds}} for 
#' details.
#'
#' @export
#' @importFrom lubridate as.duration new_interval
#' @importFrom stringr str_extract
#' @importFrom SDMTools ConnCompLabel
#' @param data_dir folder where input images are located, with filenames as 
#' output by the \code{\link{team_preprocess_landsat}} function
#' @param wrspath World Reference System (WRS) path
#' @param wrsrow World Reference System (WRS) row
#' @param start_date start date of period from which images will be chosen to 
#' fill cloudy areas in the base image (as \code{Date} object)
#' @param end_date end date of period from which images will be chosen to fill 
#' cloudy areas in the the base image (as \code{Date} object)
#' @param output_path the path to use for the output
#' @param base_date ideal date for base image (base image will be chosen as the 
#' image among the available images that is closest to this date). If NULL, 
#' then the base image will be the image with the lowest cloud cover.
#' @param fast if \code{TRUE}, use the CLOUD_REMOVE_FAST.pro script. If 
#' \code{FALSE}, use the CLOUD_REMOVE.pro script.
#' @param n_cpus the number of CPUs to use for processes that can run in 
#' parallel
#' @param overwrite whether to overwrite existing files (otherwise an error 
#' will be raised)
#' @param notify notifier to use (defaults to \code{print} function). See the 
#' \code{notifyR} package for one way of sending notifications from R. The 
#' \code{notify} function should accept a string as the only argument.
team_cloud_fill <- function(data_dir, wrspath, wrsrow, start_date, end_date, 
                            output_path, base_date=NULL, fast=FALSE, n_cpus=1, 
                            overwrite=FALSE, notify=print) {
    if (!file_test("-d", output_path)) {
        stop(paste(output_path, "does not exist"))
    }

    timer <- Track_time(notify)
    timer <- start_timer(timer, label='Cloud fill')

    if (n_cpus > 1) beginCluster(n_cpus)

    wrspath <- sprintf('%03i', wrspath)
    wrsrow <- sprintf('%03i', wrsrow)

    # Find image folders based on start and end dates
    date_regex <- '[12][8901][0-9]{2}-[123]?[0-9]{2}'
    img_dirs <- dir(data_dir,
                    pattern=paste0('^', wrspath, '-', wrsrow, '_', date_regex, 
                                   '_((LT[45])|(LE[78]))'))

    img_dates <- str_extract(img_dirs, date_regex)
    img_dates <- as.Date(img_dates, '%Y-%j')

    img_dirs <- img_dirs[(img_dates >= start_date) &
                             (img_dates < end_date)]
    img_dates <- img_dates[(img_dates >= start_date) &
                           (img_dates < end_date)]

    # Run QA stats
    masks <- list()
    imgs <- list()
    for (img_dir in img_dirs) {
        masks_file <- dir(file.path(data_dir, img_dir), pattern='masks.envi$')
        this_mask <- raster(file.path(data_dir, img_dir, masks_file), band=2)
        masks <- c(masks, this_mask)
        img_file <- dir(file.path(data_dir, img_dir), pattern='tc.envi$')
        this_img <- stack(file.path(data_dir, img_dir, img_file))
        imgs <- c(imgs, stack(this_img))
    }
    freq_table <- freq(stack(masks), useNA='no', merge=TRUE)
    # Convert frequency table to fractions
    freq_table[-1] <- freq_table[-1] / colSums(freq_table[-1], na.rm=TRUE)

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
    }

    # Convert masks to binary indicating: 0 = other; 1 = cloud or shadow
    #
    #   fmask_band key:
    #       0 = clear
    #       1 = water
    #       2 = cloud_shadow
    #       3 = snow
    #       4 = cloud
    #       255 = fill value
    for (n in 1:length(masks)) {
        masks[n] <- (masks[[n]] == 2) | (masks[[n]] == 4)
    }
    base_img <- imgs[[base_img_index]]
    imgs <- imgs[-base_img_index]
    base_mask <- masks[[base_img_index]]
    masks <- masks[-base_img_index]

    # Calculate a raster indicating the pixels in each potential fill image 
    # that are available for filling pixels of base_img that are missing due to 
    # cloud contamination. Areas coded 1 are missing due to cloud or shadow in 
    # the base image and are available in the merge image.
    fill_areas <- list()
    for (mask_img in masks) {
        fill_areas <- c(fill_areas, list(base_mask == 1 & mask_img == 0))
    }
    fill_areas_freq <- freq(stack(fill_areas), useNA='no', merge=TRUE)

    # Select the fill image with the maximum number of available pixels 
    # (counting only pixels in the fill image that are not ALSO clouded in the 
    # fill image)
    avail_fill_row <- which(fill_areas_freq$value == 1)
    fill_img_index <- which(fill_areas_freq[avail_fill_row, -1] == 
                            max(fill_areas_freq[avail_fill_row, -1]))
    fill_img <- imgs[[fill_img_index]]
    imgs <- imgs[-fill_img_index]
    fill_mask <- fill_areas[[fill_img_index]]
    masks <- masks[-base_img_index]

    # Add numbered IDs to the cloud patches using ConnCompLabel from the 
    # SDMTools package. Xiaolin Zhu's code requires these ID numbers.
    coded_cloud_mask <- ConnCompLabel(fill_mask)
    #coded_cloud_mask_patchstats <- PatchStat(coded_cloud_mask)

    cl <- options('rasterClusterObject')[[1]]
    inparallel <- FALSE
    if ((!is.null(cl)) && (nlayers(base_img) > 1)) {
        if (!require(foreach)) {
            warning('Cluster object found, but "foreach" is required to run topographic correction in parallel. Running sequentially.')
        } else if (!require(doSNOW)) {
            warning('Cluster object found, but "doSNOW" is required to run topographic correction in parallel. Running sequentially.')
        } else {
            inparallel <- TRUE
            registerDoSNOW(cl)
        }
    }

    #base_img <- dropLayer(base_img, c(4, 5, 6))
    #fill_img <- dropLayer(fill_img, c(4, 5, 6))

    if (inparallel) {
        uncorr_layer=NULL
        filled <- foreach(base_band=unstack(base_img), 
                          fill_band=unstack(fill_img), 
                          coded_cloud_mask=coded_cloud_mask, 
                          .combine='addLayer', .multicombine=TRUE, 
                          .init=raster(), .packages=c('raster', 'teamlucc', 
                                                      'rgdal')) %dopar% {
            rasterOptions(format='ENVI')
            # Write in-memory rasters to file to hand off to IDL.
            base_band <- writeRaster(base_band, rasterTmpFile(), 
                                     datatype='INT2S')
            fill_band <- writeRaster(fill_band, rasterTmpFile(), 
                                     datatype='INT2S')
            coded_cloud_mask <- writeRaster(coded_cloud_mask, rasterTmpFile(), 
                                            datatype='INT2S')
            fill_clouds(filename(base_band), filename(fill_band), 
                        filename(coded_cloud_mask), fast=fast, num_class=3, 
                        DN_min=-300, DN_max=17000, patch_long=2000)
        }
    } else {
        # Temporarily reset default raster output format to ENVI format, for 
        # compatibility with IDL.
        def_format <- rasterOptions()$format
        rasterOptions(format='ENVI')
        base_img <- writeRaster(base_img, rasterTmpFile(), 
                                datatype='INT2S')
        fill_img <- writeRaster(fill_img, rasterTmpFile(), 
                                 datatype='INT2S')
        coded_cloud_mask <- writeRaster(coded_cloud_mask, rasterTmpFile(), 
                                        datatype='INT2S')
        filled <- fill_clouds(filename(base_img), filename(fill_img), 
                              filename(coded_cloud_mask), fast=fast, 
                              num_class=3, DN_min=-300, DN_max=17000, 
                              patch_long=2000)
        rasterOptions(format=def_format)
    }

    timer <- stop_timer(timer, label='Cloud fill')


    if (n_cpus > 1) endCluster()

    return(filled)
}
