\name{track_time}
\alias{track_time}
\title{Print status messages tracking run time for an R script}
\usage{
track_time(label = "", action = "print")
}
\arguments{
  \item{label}{an optional label used to maintain multiple
  tracking timers}

  \item{action}{a string indicating whether to \code{start}
  the timer or \code{print} (default) the time that has
  passed since the timer was last started.}
}
\description{
Print status messages tracking run time for an R script
}
\examples{
track_time(action="start")
track_time()

track_time(label='process1', action="start")
track_time(label='process1')

track_time()
}

