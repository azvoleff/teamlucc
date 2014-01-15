\name{startTimer}
\alias{startTimer}
\title{Start a tracking timer}
\usage{
startTimer(x, label)
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

timer <- stopTimer('test')
}
\seealso{
\code{\link{stopTimer}}
}

