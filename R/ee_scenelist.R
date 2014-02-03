#' Save scenelist from EarthExplorer metadata for upload to ESPA
#'
#' @export
#' @importFrom lubridate as.interval new_period %within%
#' @param x a \code{data.frame} with a list of Landsat scenes as output from 
#' the save metadata function on http://earthexplorer.usgs.gov
#' @param start_date starting date as a \code{Date} object
#' @param out_file filename for output text file for later upload to ESPA
#' @param num_months number of months to include
#' @param min_clear the minimum percent clear to plot (calculated as 1 - 
#' percent cloud cover). Images with less than \code{min_clear} fraction of the 
#' image area free of clouds will be ignored.
#' @param exclude a list of sensors to exclude (for example, set 
#' \code{exclude=c('LE7', 'LT4')} to exclude Landsat 7 ETM+ and Landsat 4 TM 
#' images.
#' @return used for side effect of producing ESPA scene list
ee_scenelist <- function(x, start_date, out_file, num_months=12,
                         min_clear=.7, exclude=list()) {
    if (!class(start_date) == 'Date') {
        stop('start_date must be a "Date" object')
    }
    sel_interval <- as.interval(new_period(months=num_months), start_date)
    x$Date.Acquired <- as.Date(as.character(x$Date.Acquired))
    x <- x[x$Date.Acquired %within% sel_interval, ]
    x$Sensor <- substr(x$Landsat.Scene.Identifier, 1, 3)
    x <- x[!(x$Sensor %in% exclude), ]
    x$Frac_Clear <- (100 - x$Cloud.Cover) / 100
    x <- x[x$Frac_Clear >= min_clear, ]
    write.table(x$Landsat.Scene.Identifier, out_file, row.names=FALSE, 
                col.names=FALSE, quote=FALSE, sep='\n')
}
