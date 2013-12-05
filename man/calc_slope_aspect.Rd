\name{calc_slope_aspect}
\alias{calc_slope_aspect}
\title{Calculates slope and aspect using same formula as in the "landsat"
package by Sarah Goslee.}
\usage{
  calc_slope_aspect(Rast, EWkernel, EWres, NSkernel, NSres,
    smoothing)
}
\arguments{
  \item{Rast}{a matrix of elevations. The slope and aspect
  will be calculated for the center pixel in the matrix.}

  \item{EWkernel}{kernel to use for East-West gradient
  calculation}

  \item{EWres}{East-West resolution}

  \item{NSkernel}{kernel to use for North-South gradient
  calculation}

  \item{NSres}{North-South resolution}

  \item{smoothing}{positive integer for smoothing. 1 means
  no smoothing.}
}
\value{
  list with two elements. The first element is the
  calculated slope, the second element is the calculated
  aspect.
}
\description{
  This function is called by the \code{\link{slopeasp_par}}
  function. It is not intended to be used directly.
}
\references{
  Sarah Goslee. Analyzing Remote Sensing Data in {R}: The
  {landsat} Package. Journal of Statistical Software, 2011,
  43:4, pg 1--25. http://www.jstatsoft.org/v43/i04/
}

