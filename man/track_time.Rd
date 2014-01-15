\name{Track_time}
\alias{Track_time}
\title{Instantiate a new Track_time object}
\usage{
Track_time(notify = print)
}
\arguments{
  \item{label}{an optional label used to maintain multiple
  tracking timers}
}
\description{
Creates a new Track_time object for use in tracking and
printing status the running time of processes in an R
script.
}
\examples{
timer <- Track_time()
print(timer)
}
\seealso{
\code{\link{startTimer}}, \code{\link{stopTimer}}
}

