# Function for retrieving frequencies from a raster frequency table, given a 
# band name and value. Handles the case of a given code not occurring in a 
# particular raster, in which case it will not show up as a row in the 
# frequency table.
get_freq <- function(band, value, freq_table) {
    band_col <- grep(band, names(freq_table))
    if (!(value %in% freq_table[, 1])) {
        # Return 0 if this value never shows up in the raster
        return(0)
    }
    frac <- freq_table[freq_table[1] == value, band_col]
    if (is.na(frac)) {
        return(0)
    } else {
        return(round(frac, 4))
    }
}

#' Calculate statistics on imagery within an AOI
#'
#' @export
#' @importFrom stringr str_extract
#' @param image_dirs list of paths to a set of Landsat CDR image files in 
#' GeoTIFF format as output by the \code{unstack_ledapscdr} function.
#' @param aoi an area of interest (AOI) to crop from each image
#' @return a \code{data.frame}
auto_QA_stats <- function(image_dirs, aoi) {
    lndsr_regex <- '^(lndsr.)?((LT4)|(LT5)|(LE7)|(LC8))[0-9]{6}[12][0-9]{6}[a-zA-Z]{3}[0-9]{2}'
    mask_bands <- c('fill_QA', 'fmask_band')

    out <- c()
    for (image_dir in image_dirs) {
        lndsr_files <- dir(image_dir, pattern=lndsr_regex)
        image_basenames <- unique(str_extract(lndsr_files,lndsr_regex))

        if (length(image_basenames) == 0) {
            stop(paste('no files found in', image_dir))
        }

        for (image_basename in image_basenames) {
            message(paste0('Processing ', image_basename, '...'))
            metadata_string <- str_extract(image_basename, 
                                           '((LT4)|(LT5)|(LE7)|(LC8))[0-9]{13}')
            sensor <- str_extract(metadata_string, '^((LT[45])|(LE7)|(LC8))')
            year <- substr(metadata_string, 10, 13)
            julian_day <- substr(metadata_string, 14, 16)
            img_path <- substr(metadata_string, 4, 6)
            img_row <- substr(metadata_string, 7, 9)

            mask_band_files <- c()
            for (mask_band in mask_bands) {
                mask_band_files <- c(mask_band_files,
                                     paste(file.path(image_dir, 
                                                     image_basename), 
                                           mask_band, sep='_'))
            }
            mask_band_files <- paste0(mask_band_files, '.tif')
            mask_stack <- stack(mask_band_files)
            names(mask_stack) <- mask_bands
            if (!missing(aoi)) {
                if (proj4string(aoi) != proj4string(mask_stack)) {
                    if (class(aoi) == 'Raster') {
                        aoi <- projectRaster(aoi, mask_stack)
                    } else {
                        aoi <- spTransform(aoi, CRS(proj4string(mask_stack)))
                    }
                }
                mask_stack <- crop(mask_stack, aoi)
                mask_stack <- mask(mask_stack, aoi)
            }

            freq_table <- freq(mask_stack, useNA='no', merge=TRUE)
            # Convert frequency table to fractions
            freq_table[-1] <- freq_table[-1] / colSums(freq_table[-1], na.rm=TRUE)

            out <- c(out, list(list(img_path,
                                    img_row,
                                    year,
                                    julian_day,
                                    sensor,
                                    get_freq('fill_QA', 0, freq_table),
                                    get_freq('fill_QA', 255, freq_table),
                                    get_freq('fmask_band', 0, freq_table),
                                    get_freq('fmask_band', 1, freq_table),
                                    get_freq('fmask_band', 2, freq_table),
                                    get_freq('fmask_band', 3, freq_table),
                                    get_freq('fmask_band', 4, freq_table),
                                    get_freq('fmask_band', 255, freq_table))))
        }
    }

    out <- data.frame(matrix(unlist(out), nrow=length(out), byrow=T))
    names(out) <- c('path', 'row', 'year', 'julian', 'sensor',
                    'fill_QA_notfill', 'fill_QA_fill', 'fmask_clear', 
                    'fmask_water', 'fmask_cloud_shadow', 'fmask_snow', 
                    'fmask_cloud', 'fmask_fill')
    return(out)
}
