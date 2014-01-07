#' Print status messages tracking run time for an R script
#'
#' @export
#' @importFrom lubridate as.duration now
#' @param label an optional label used to maintain multiple tracking timers
#' @param action a string indicating whether to \code{start} the timer or 
#' \code{print} (default) the time that has passed since the timer was last 
#' started.
#' @examples
#' track_time(action="start")
#' track_time()
#'
#' track_time(label='process1', action="start")
#' track_time(label='process1')
#'
#' track_time()
track_time <- function(label='', action='print') {
    if (action == 'start') {
        assign(paste('track_time_starttime', label, sep='_'), now(), 
               inherits=TRUE)
        if (label != '') {
            return(paste0('track_time: Tracking timer "',
                          label,'" started at ', now(), '.'))
        } else {
            return(paste0('track_time: Tracking timer started at ',
                          now(), '.'))
        }
    } else if (action == 'print') {
        startTimeVar <- mget(paste('track_time_starttime', label, sep='_'), 
                             ifnotfound=list(NULL), inherits=TRUE)
        if (!is.null(startTimeVar[[1]])) {
            if (label != '') {
                return(paste0('track_time (', label, '): ',
                              round(as.duration(now() - startTimeVar[[1]]), 2)))
            } else {
                return(paste('track_time:',
                             round(as.duration(now() - startTimeVar[[1]]), 2)))
            }
        } else {
            if (label == '') {
                stop(paste0('No unlabeled timer found. Must call track_time(action="start") before calling track_time()'))
            } else {
                stop(paste0('No timer with label "', label, '" found. Must call track_time(action="start") before calling track_time()'))
            }
        }
    } else {
        stop(paste('Unrecognized action "', action, '"', sep=''))
    }
}
