\name{stopTimer}
\alias{stopTimer}
\title{Stop a tracking timer}
\usage{
stopTimer(x, label = "Default")
}
\arguments{
  \item{label}{an optional label used to maintain multiple
  tracking timers}
}
\description{
The \code{label} is optional. If not supplied, the default
timer, labelled 'Default' will be used.
}
\examples{
timer <- Track_time()
print(timer)

timer <- startTimer('test')

print(timer)

print(timer, 'test')
}
\seealso{
\code{\link{startTimer}}
}

