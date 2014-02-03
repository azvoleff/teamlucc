#' Plot EarthExplorer metadata file
#'
#' @export
#' @import ggplot2
#' @importFrom plyr ddply .
#' @param x a \code{data.frame} with a list of Landsat scenes as output from 
#' the save metadata function on http://earthexplorer.usgs.gov
#' @param start_year an integer indicating the first year to plot
#' @param end_year an integer indicating the last year to plot
#' @param min_clear the minimum percent clear to plot (calculated as 1 - 
#' percent cloud cover). Images with less than \code{min_clear} fraction of the 
#' image area free of clouds will be ignored.
#' @param exclude a list of sensors to exclude (for example, set 
#' \code{exclude=c('LE7', 'LT4')} to exclude Landsat 7 ETM+ and Landsat 4 TM 
#' images.
#' @return used for side effect of producing a plot
plot_ee <- function(x, start_year, end_year, min_clear=.7, exclude=list()) {
    x$Sensor <- substr(x$Landsat.Scene.Identifier, 1, 3)
    x <- x[!(x$Sensor %in% exclude), ]
    x$Sensor <- factor(x$Sensor)
    x$Date.Acquired <- as.Date(as.character(x$Date.Acquired))
    x$Year <- as.numeric(format(x$Date.Acquired, '%Y'))
    x$Month <- as.numeric(format(x$Date.Acquired, '%m')) - .5
    x$MonthFactor <- factor(format(x$Date.Acquired, '%m'))
    x$Path_Row <- factor(paste(x$WRS.Path, x$WRS.Row, sep='/'))
    x$YearMonth <- paste(x$Year, x$MonthFactor, sep='/')
    x$Frac_Clear <- (100 - x$Cloud.Cover) / 100
    x <- x[order(x$WRS.Path, x$WRS.Row), ]
    x <- x[x$Frac_Clear >= min_clear, ]
    if (!missing(start_year)) {
        x <- x[x$Year >= start_year, ]
    }
    if (!missing(end_year)) {
        x <- x[x$Year <= end_year, ]
    }
    if (nrow(x) == 0) {
        stop('no data to plot - try different start/end years')
    }
    YearMonth=Month=Cum_Month=Path_Row=Sensor=Frac_Clear=NULL # Keep R CMD CHECK happy
    x <- ddply(x, .(YearMonth), transform,
               Cum_Month=cumsum(rep(1, length(Month))))
    p <- ggplot(x, aes(xmin=Month,
                       xmax=Month + 1, 
                       ymin=Cum_Month - 1, 
                       ymax=Cum_Month,
                       colour=Sensor,
                       fill=Path_Row,
                       alpha=Frac_Clear)) +
        geom_rect() + facet_grid(Year ~ ., scales='free_y', space='free_y') +
        xlab('Month') +
        scale_x_continuous(breaks=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12),
                           labels=c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 
                                    'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec')) +
        theme(panel.grid.minor.y=element_blank(), 
              panel.grid.major.y=element_blank(),
              panel.grid.major.x=element_blank()) +
        theme(axis.ticks.y=element_blank(),
              axis.text.y=element_blank()) +
        scale_colour_brewer(type='qual', palette='Set1', drop=FALSE, name='Sensor') +
        scale_fill_brewer(type='qual', palette='Set2', drop=FALSE, name='Path/Row') +
        scale_alpha(name='Fraction Clear')
    return(p)
}
