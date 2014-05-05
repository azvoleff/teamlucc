#' Save scenelist from EarthExplorer metadata for upload to ESPA
#'
#' @export
#' @importFrom lubridate new_interval %within%
#' @param x a \code{data.frame} with a list of Landsat scenes as output from 
#' the save metadata function on http://earthexplorer.usgs.gov
#' @param start_date starting date as a \code{Date} object
#' @param end_date end date as a \code{Date} object
#' @param out_file filename for output text file for later upload to ESPA
#' @param min_clear the minimum percent clear to plot (calculated as 1 - 
#' percent cloud cover). Images with less than \code{min_clear} fraction of the 
#' image area free of clouds will be ignored.
#' @param exclude a list of sensors to exclude (for example, set 
#' \code{exclude=c('LE7', 'LT4')} to exclude Landsat 7 ETM+ and Landsat 4 TM 
#' images.
#' @return used for side effect of producing ESPA scene list
espa_scenelist <- function(x, start_date, end_date, out_file, min_clear=.7, 
                         exclude=list()) {
    if (!class(start_date) == 'Date') {
        stop('start_date must be a "Date" object')
    }
    if (!class(end_date) == 'Date') {
        stop('end_date must be a "Date" object')
    }
    if ((!missing(start_date) && missing(end_date)) ||
        (missing(start_date) && !missing(end_date))) {
        stop('both start_date and end_date must be provided')
    } else if (!missing(start_date) && !missing(end_date)) {
        sel_interval <- new_interval(start_date, end_date)
        x <- x[x$Date.Acquired %within% sel_interval, ]
    }
    if (nrow(x) == 0) {
        stop('no data to download - try different start/end dates')
    }
    x <- x[!(x$Sensor %in% exclude), ]
    x <- x[x$Frac_Clear >= min_clear, ]
    write.table(x$Landsat.Scene.Identifier, out_file, row.names=FALSE, 
                col.names=FALSE, quote=FALSE, sep='\n')
}
