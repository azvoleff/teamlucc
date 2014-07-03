#' A class for tracking running time of individual sections of an R script
#' @slot timers a \code{data.frame} tracking timer names and start times
#' @slot notify function to use for outputting timers (defaults to 
#' \code{\link{print}}
#' @import methods
#' @importFrom lubridate now
#' @export Track_time
#' @name Track_time-class
setClass('Track_time', slots=c(timers='data.frame', notify="function"),
    prototype=list(timers=data.frame(label='Default', starttime=now()), 
                   notify=print)
)

#' Instantiate a new Track_time object
#'
#' Creates a new Track_time object for use in tracking and printing status the 
#' running time of processes in an R script.
#'
#' @export Track_time
#' @import methods
#' @importFrom lubridate now
#' @param notify a function to handle the string output from Track_time.  This 
#' function should accept a string as an argument. Default is the
#' \code{\link{print}} function.
#' @return Track_time object
#' @seealso \code{\link{start_timer}}, \code{\link{stop_timer}}
#' @examples
#' timer <- Track_time()
#' print(timer)
Track_time <- function(notify=print) {
    return(new('Track_time',
               timers=data.frame(label='Default', starttime=now()),
               notify=notify))
}

#' Print a Track_time object
#'
#' @export
#' @importFrom lubridate now as.duration
#' @import methods
#' @param x a Track_time object
#' @param label (optional) selects a specific tracking timer to print
#' @param ... ignored
#' @method print Track_time
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
.start_timer <- function(x, label) {
    if (!missing(label)) {
        if (label %in% x@timers$label) {
            stop(paste0('"', label, '"', ' timer already defined'))
        }
        x@timers <- rbind(x@timers, data.frame(label=label, starttime=now()))
        x@notify(paste0(x@timers$starttime[x@timers$label == label], ': started "', label, '"'))
    } else {
        x@timers$starttime[x@timers$label == 'Default'] <- now()
        x@notify(paste0(x@timers$starttime[x@timers$label == "Default"], ': started'))
    }
    return(x)
}

#' Start a tracking timer
#'
#' The \code{label} is optional. If not supplied the default timer (labelled 
#' "Default") will be used.
#'
#' @export start_timer
#' @param x a \code{Track_time} object
#' @param label an optional label used to maintain multiple tracking timers
#' @return Track_time object
#' @seealso \code{\link{stop_timer}}
#' @examples
#' timer <- Track_time()
#' print(timer)
#'
#' timer <- start_timer(timer, 'test')
#'
#' print(timer, 'test')
#'
#' timer <- stop_timer(timer, 'test')
#'
#' print(timer)
setGeneric("start_timer", function(x, label) {
    standardGeneric("start_timer")
})

#' @rdname start_timer
#' @aliases start_timer,Track_time-method
setMethod("start_timer", signature(x="Track_time"),
    function(x) .start_timer(x)
)

#' @rdname start_timer
#' @aliases start_timer,Track_time,character-method
setMethod("start_timer", signature(x="Track_time", label="character"),
    function(x, label) .start_timer(x, label)
)

#' @importFrom lubridate now as.duration
.stop_timer <- function(x, label='Default') {
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
    x@notify(paste0(now(), ': finished "', label, '" (', round(elapsed, 3),' elapsed)'))
    return(x)
}

#' Stop a tracking timer
#'
#' The \code{label} is optional. If not supplied, the default timer, labelled 
#' 'Default' will be used.
#'
#' @export stop_timer
#' @param x a \code{Track_time} object
#' @param label an optional label used to maintain multiple tracking timers
#' @seealso \code{\link{start_timer}}
#' @return Track_time object
#' @examples
#' timer <- Track_time()
#' print(timer)
#'
#' timer <- start_timer(timer, 'test')
#'
#' print(timer, 'test')
#'
#' timer <- stop_timer(timer, 'test')
#'
#' print(timer)
setGeneric("stop_timer", function(x, label='Default') {
    standardGeneric("stop_timer")
})

#' @rdname stop_timer
#' @aliases stop_timer,Track_time-method
setMethod("stop_timer", signature(x="Track_time"),
    function(x) .stop_timer(x)
)

#' @rdname stop_timer
#' @aliases stop_timer,Track_time,character-method
setMethod("stop_timer", signature(x="Track_time", label="character"),
    function(x, label) .stop_timer(x, label)
)
