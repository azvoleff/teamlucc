\name{start_timer}
\alias{start_timer}
\title{Start a tracking timer}
\usage{
start_timer(x, label)
}
\arguments{
  \item{x}{a \code{Track_time} object}

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

timer <- start_timer(timer, 'test')

print(timer, 'test')

timer <- stop_timer(timer, 'test')

print(timer)
}
\seealso{
\code{\link{stop_timer}}
}

