\name{slopeasp_par}
\alias{slopeasp_par}
\title{Slope and aspect calculation in parallel for DEM \code{RasterLayer}}
\usage{
slopeasp_par(DEM, EWkernel, NSkernel, smoothing = 1, filename = NULL, ...)
}
\arguments{
  \item{DEM}{DEM as \code{RasterLayer}}

  \item{EWkernel}{kernel to use for East-West gradient
  calculation}

  \item{NSkernel}{kernel to use for North-South gradient
  calculation}

  \item{smoothing}{positive integer for smoothing. 1 means
  no smoothing.}

  \item{filename}{output file for 2-band slope and aspect
  layer stack (optional)}

  \item{...}{additional parameters to pass to rasterEngine}
}
\value{
RasterBrick with two layers: 'slope' and 'aspect'
}
\description{
This code calculate the slope and aspect for a DEM provided
as a \code{RasterLayer}. The computations will be done in
parallel if \code{sfQuickInit()} is called before
\code{slopeasp_par}.
}
\details{
The code for calculating the slope and aspect is taken
directly from the \code{landsat} package by Sarah Goslee.
See the help page for \code{slopeasp} in the \code{landsat}
package for details on the parameters.
}
\examples{
#TODO: Write examples
}
\author{
Sarah Goslee and Alex Zvoleff
}
\references{
Sarah Goslee. Analyzing Remote Sensing Data in {R}: The
{landsat} Package. Journal of Statistical Software, 2011,
43:4, pg 1--25. http://www.jstatsoft.org/v43/i04/
}

