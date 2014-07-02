#' Plot EarthExplorer scene list
#'
#' This function can produce two different types of plots from a USGS 
#' EarthExplorer Landsat CDR Surface Reflectance scene list.
#'
#' @export
#' @import ggplot2
#' @importFrom dplyr group_by summarize
#' @importFrom lubridate new_interval %within%
#' @param x a \code{data.frame} with a list of Landsat scenes as output by
#' \code{\link{ee_read}}
#' @param start_date starting date as a \code{Date} object
#' @param end_date end date as a \code{Date} object
#' @param min_clear the minimum percent clear to plot (calculated as 1 - 
#' percent cloud cover). Images with less than \code{min_clear} fraction of the 
#' image area free of clouds will be ignored.
#' @param exclude a list of sensors to exclude (for example, set 
#' \code{exclude=c('LE7', 'LT4')} to exclude Landsat 7 ETM+ and Landsat 4 TM 
#' images.
#' @param normalize if \code{TRUE}, plot as a normalized line plot
#' @param title title for plot (or \code{NULL} for no title)
#' @return used for side effect of producing a plot
ee_plot <- function(x, start_date, end_date, min_clear=.7, exclude=list(), 
                    normalize=FALSE, title=NULL) {
    if (!class(start_date) == 'Date') {
        stop('start_date must be a "Date" object')
    }
    if (!class(end_date) == 'Date') {
        stop('end_date must be a "Date" object')
    }
    x <- x[!(x$Sensor %in% exclude), ]
    x$Sensor <- factor(x$Sensor)

    x <- x[order(x$WRS.Path, x$WRS.Row), ]
    x <- x[x$Frac_Clear >= min_clear, ]

    if ((!missing(start_date) && missing(end_date)) ||
        (missing(start_date) && !missing(end_date))) {
        stop('both start_date and end_date must be provided')
    } else if (!missing(start_date) && !missing(end_date)) {
        sel_interval <- new_interval(start_date, end_date)
        x <- x[x$Date.Acquired %within% sel_interval, ]
    }
    if (nrow(x) == 0) {
        stop('no data to plot - try different start/end dates')
    }

    if (!normalize) {
        YearMonth=Month=Cum_Month=Path_Row=Sensor=Frac_Clear=NULL # Keep R CMD CHECK happy
        x <- transform(group_by(x, YearMonth),
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
    } else {
        # Keep R CMD CHECK happy:
        YearMonth=Path_Row=Year=Month=Max_Frac_Clear=Frac_Clear=Sum_Max_Frac_Clear=NULL
        Frac_Clear_Stats <- summarize(group_by(x, YearMonth, Path_Row),
                                      Year=Year[1], Month=Month[1],
                                      Max_Frac_Clear=max(Frac_Clear))
        Frac_Clear_Stats <- summarize(group_by(Frac_Clear_Stats, YearMonth), 
                                      Year=Year[1], Month=Month[1],
                                      Sum_Max_Frac_Clear=sum(Max_Frac_Clear))
        p <- ggplot(Frac_Clear_Stats, aes(xmin=Month + .05,
                                          xmax=Month + 1-.05, 
                                          ymin=0, 
                                          ymax=Sum_Max_Frac_Clear)) +
            geom_rect() + facet_grid(Year ~ .) +
            xlab('Month') + ylab('Total Fraction Clear') +
            scale_x_continuous(breaks=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12),
                               labels=c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 
                                        'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec')) +
            theme(panel.grid.minor.y=element_blank(), 
                  panel.grid.major.y=element_blank(),
                  panel.grid.major.x=element_blank()) +
            theme(axis.ticks.y=element_blank(),
                  axis.text.y=element_blank()) +
            geom_hline(yintercept=seq(1,length(unique(x$Path_Row))), colour='white', linetype='dashed')
    }

    if (!is.null(title)) {
        p <- p + ggtitle(title)
    }

    return(p)
}
