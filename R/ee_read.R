#' Read EarthExplorer CSV format scene list
#'
#' This function reads in a CSV file of Landsat CDR Surface Reflectance images 
#' as output from USGS EarthExplorer. param x a \code{data.frame} with a list 
#' of Landsat scenes as output from the save metadata function on 
#'
#' @export
#' @param x path to a CSV file with a list of Landsat scenes as output from the 
#' save metadata function on http://earthexplorer.usgs.gov
#' @return x a \code{data.frame} with a list of Landsat scenes and their 
#' associated metadata
ee_read <- function(x) {
    scenes <- read.csv(x, stringsAsFactors=FALSE, quote="", 
                          na.strings=c('NA', ' '))

    scenes$Sensor <- substr(scenes$Landsat.Scene.Identifier, 1, 3)
    scenes$Sensor <- factor(scenes$Sensor)

    # Dates are formatted as either: 1999/12/31 or 12/31/1999
    yr_first <- grepl('^[0-9]{4}/', scenes$Date.Acquired)
    yr_last <- grepl('/[0-9]{4}$', scenes$Date.Acquired)
    if ((sum(yr_first) + sum(yr_last)) < nrow(scenes)) {
        stop('unrecognized date format in Date.Acquired column')
    }
    acq_date <- as.Date(scenes$Date.Acquired)
    acq_date[yr_first] <- as.Date(scenes$Date.Acquired[yr_first], '%Y/%m/%d')
    acq_date[yr_last] <- as.Date(scenes$Date.Acquired[yr_last], '%m/%d/%Y')
    scenes$Date.Acquired <- acq_date

    scenes$Year <- as.numeric(format(scenes$Date.Acquired, '%Y'))
    scenes$Month <- as.numeric(format(scenes$Date.Acquired, '%m')) - .5
    scenes$MonthFactor <- factor(format(scenes$Date.Acquired, '%m'))
    scenes$Path_Row <- factor(paste(scenes$WRS.Path, scenes$WRS.Row, sep='/'))
    scenes$YearMonth <- paste(scenes$Year, scenes$MonthFactor, sep='/')
    scenes$Frac_Clear <- (100 - scenes$Cloud.Cover) / 100
    scenes <- scenes[order(scenes$WRS.Path, scenes$WRS.Row), ]

    if (nrow(scenes) == 0) {
        stop(paste0('no data found in', x,
                    ' - is scenes an EarthExplorer CSV export?'))
    }

    # Drop ENGINEERING, TEST, EXCHANGE, and VALIDATION data. See the Landsat 
    # data dictionary at: https://lta.cr.usgs.gov/landsat_dictionary.html
    scenes <- scenes[scenes$Data.Category == 'NOMINAL', ]

    return(scenes)
}
