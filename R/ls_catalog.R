#' Catalog a folder of Landsat images
#'
#' This function is used to produce a \code{data.frame} of Landsat images 
#' stored locally after download from the USGS. The images should be in a 
#' series of subfolders named following the naming scheme of 
#' \code{\link{espa_extract}}.
#'
#' @export
#' @importFrom stringr str_extract
#' @param in_folder path to a folder of Landsat surface reflectance images (for 
#' example, as extracted by the \code{espa_extract} function).
#' @return a \code{data.frame} with a list of the Landsat images found within 
#' in_folder
ls_catalog <- function(in_folder) {
    if (!file_test('-d', in_folder)) {
        stop(paste(in_folder, 'does not exist'))
    }
    out <- c()
    for (outer_item in dir(in_folder)) {
        outer_item_full <- file.path(in_folder, outer_item)
        # First cycle through the yearly folders
        if (!file_test('-d', outer_item_full)) {
            next
        }
        for (inner_item in dir(outer_item_full)) {
            inner_item_full <- file.path(outer_item_full, inner_item)
            # Check to ensure inner item is a Landsat HDF file
            if (!file_test('-f', inner_item_full) ||
                !grepl('^(lndsr.)?((LT4)|(LT5)|(LE7)|(LC8))[0-9]{13}[A-Z]{3}[0-9]{2}.hdf$', 
                       inner_item)) {
                next
            }
            metadata_string <- str_extract(inner_item, 
                                           '((LT4)|(LT5)|(LE7)|(LC8))[0-9]{13}')
            if (grepl('^LT4', metadata_string)) {
                sensor <- 'LT4'
            } else if (grepl('^LT5', metadata_string)) {
                sensor <- 'LT5'
            } else if (grepl('^LE7', metadata_string)) {
                sensor <- 'LE7'
            } else if (grepl('^LC8', metadata_string)) {
                sensor <- 'LC8'
            } else {
                message(paste('Skipping', inner_item,
                              '- cannot determine sensor from filename.'))
                next
            }
            year <- substr(metadata_string, 10, 13)
            julian_day <- substr(metadata_string, 14, 16)
            img_path <- substr(metadata_string, 4, 6)
            img_row <- substr(metadata_string, 7, 9)
            img_date <- as.Date(paste0(year, julian_day), '%Y%j')
            out <- c(out, list(list(img_path,
                                    img_row,
                                    year,
                                    julian_day,
                                    sensor,
                                    format(img_date, '%m'), 
                                    format(img_date, '%d'),
                                    inner_item)))

        }
    }
    if (is.null(out)) {
        stop('no images found')
    }
    out <- data.frame(matrix(unlist(out), nrow=length(out), byrow=T))
    names(out) <- c('path', 'row', 'year', 'julian', 'sensor',
                    'month', 'day', 'filename')
    out <- out[order(out$path, out$row, out$year, out$julian, out$sensor), ]
    return(out)
}
