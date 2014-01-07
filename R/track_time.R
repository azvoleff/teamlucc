#' Print status messages tracking run time for an R script
#'
#' @export
#' @importFrom lubridate as.duration now
#' @param label an optional label used to maintain multiple tracking timers
#' @param action a string indicating whether to \code{start} the timer or 
#' \code{print} (default) the time that has passed since the timer was last 
#' started.
#' @examples
#' #TODO: Add examples
track_time <- function(label='', action='print') {
    if (action == 'start') {
        assign(paste('trackTime_starttime', label, sep='_'), now(), 
               inherits=TRUE)
        if (label != '') {
            return(paste0('trackTime: Tracking timer "', label,'" started at ', now(), '.'))
        } else {
            return(paste0('trackTime: Tracking timer started at ', now(), '.'))
        }
    } else if (action == 'print') {
        startTimeVar <- mget(paste('trackTime_starttime', label, sep='_'), 
                             ifnotfound=list(NULL), inherits=TRUE)
        if (!is.null(startTimeVar[[1]])) {
            if (label != '') {
                return(paste0('trackTime (', label, '): ', round(as.duration(now() - startTimeVar[[1]]), 2)))
            } else {
                return(paste('trackTime:', round(as.duration(now() - startTimeVar[[1]]), 2)))
            }
        } else {
            stop(paste0('No timer with label "', label, '" found. Must call trackTime(action="start") before calling trackTime()'))
        }
    } else {
        stop(paste('Unrecognized action "', action, '"', sep=''))
    }
}
