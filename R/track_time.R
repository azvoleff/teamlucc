#' Class to represent a tracking timer
#'
#' @import methods
#' @importFrom lubridate now
#' @export
#' @name track_time-class
setClass('Track_time', slots=c(timers='data.frame', notify="function"),
    prototype=list(timers=data.frame(label='Default', starttime=now()), 
                   notify=print)
)

#' Instantiate a new Track_time object
#'
#' Creates a new Track_time object for use in tracking and printing status the 
#' running time of processes in an R script.
#'
#' @export
#' @import methods
#' @importFrom lubridate now
#' @param label an optional label used to maintain multiple tracking timers
#' @seealso \code{\link{startTimer}}, \code{\link{stopTimer}}
#' @examples
#' timer <- Track_time()
#' print(timer)
Track_time <- function(notify=print) {
    return(new('Track_time',
               timers=data.frame(label='Default', starttime=now()),
               notify=notify))
}

#' @import methods
#' @importFrom lubridate now as.duration
#' @S3method print Track_time
print.Track_time <- function(x, label, ...) {
    timers <- x@timers
    if (!missing(label)) {
        if (!(label %in% timers$label)) {
            stop(paste0('"', label, '"', ' timer not defined'))
        } else {
            timers <- timers[timers$label == label, ] 
        }
    }
    for (n in 1:nrow(timers)) {
       x@notify(paste(timers$label[n], 'timer:',
                round(as.duration(now() - timers$starttime[n]), 3),
                'elapsed'))
    }
}

setMethod("show", signature(object="Track_time"), function(object) print(object))

#' @importFrom lubridate now as.duration
.startTimer <- function(x, label) {
    if (!missing(label)) {
        if (label %in% x@timers$label) {
            stop(paste0('"', label, '"', ' timer already defined'))
        }
        x@timers <- rbind(x@timers, data.frame(label=label, starttime=now()))
        x@notify(paste(label, 'timer started at:', 
                        x@timers$starttime[x@timers$label == label]))
    } else {
        x@timers$starttime[x@timers$label == 'Default'] <- now()
        x@notify(paste('Default timer started at:', 
                       x@timers$starttime[x@timers$label == 'Default']))
    }
    return(x)
}

#' Start a tracking timer
#'
#' The \code{label} is optional. If not supplied, the default timer, labelled 
#' 'Default' will be used.
#'
#' @export
#' @param label an optional label used to maintain multiple tracking timers
#' @seealso \code{\link{stopTimer}}
#' @examples
#' timer <- Track_time()
#' print(timer)
#'
#' timer <- startTimer('test')
#'
#' print(timer)
#'
#' print(timer, 'test')
#'
#' timer <- stopTimer('test')
setGeneric("startTimer", function(x, label) {
    standardGeneric("startTimer")
})

setMethod("startTimer", signature(x="Track_time"),
    function(x) .startTimer(x)
)

setMethod("startTimer", signature(x="Track_time"),
    function(x, label) .startTimer(x, label)
)

#' @importFrom lubridate now as.duration
.stopTimer <- function(x, label='Default') {
    if (!(label %in% x@timers$label)) {
        stop(paste0('"', label, '"', ' timer not defined'))
    }
    elapsed <- as.duration(now() - x@timers$starttime[x@timers$label == label])
    if (label == 'Default') {
        # Never delete the default timer. Only reset it.
        x@timers$starttime[x@timers$label == 'Default'] <- now()
    } else {
        x@timers <- x@timers[x@timers$label != label, ] 
    }
    x@notify(paste0(label, ' timer finished at: ', now(), ' (', round(elapsed, 3),' elapsed)'))
    return(x)
}

#' Stop a tracking timer
#'
#' The \code{label} is optional. If not supplied, the default timer, labelled 
#' 'Default' will be used.
#'
#' @export
#' @param label an optional label used to maintain multiple tracking timers
#' @seealso \code{\link{startTimer}}
#' @examples
#' timer <- Track_time()
#' print(timer)
#'
#' timer <- startTimer('test')
#'
#' print(timer)
#'
#' print(timer, 'test')
#'
setGeneric("stopTimer", function(x, label='Default') {
    standardGeneric("stopTimer")
})

setMethod("stopTimer", signature(x="Track_time"),
    function(x) .stopTimer(x)
)

setMethod("stopTimer", signature(x="Track_time"),
    function(x, label) .stopTimer(x, label)
)
