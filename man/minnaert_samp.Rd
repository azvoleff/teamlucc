\name{minnaert_samp}
\alias{minnaert_samp}
\title{Topographic correction for satellite imagery using Minnaert method}
\usage{
  minnaert_samp(x, slope, aspect, sunelev, sunazimuth,
    na.value = NA, GRASS.aspect = FALSE,
    IL.epsilon = 1e-06,
    slopeclass = c(1, 5, 10, 15, 20, 25, 30, 45),
    coverclass = NULL, sampleindices = NULL)
}
\arguments{
  \item{x}{image as a \code{RasterLayer}}

  \item{slope}{the slope as a \code{RasterLayer}}

  \item{aspect}{the aspect as a \code{RasterLayer}}

  \item{sunelev}{sun elevation in degrees}

  \item{sunazimuth}{sun azimuth in degrees}

  \item{na.value}{the value used to code no data values}

  \item{GRASS.aspect}{is aspect defined according to GRASS
  convetion}

  \item{IL.epsilon}{a small amount to add to calculated
  illumination values that are equal to zero to avoid
  division by zero resulting in Inf values}

  \item{slopeclass}{the slope classes to calculate k for}

  \item{coverclass}{used to calculate k for specific cover
  class (optional)}

  \item{sampleindices}{(optional) row-major indices of
  sample pixels to use in the calculation of k values for
  the Minnaert correction. See \code{\link{gridsample}}.}
}
\value{
  \code{RasterLayer} with topographically corrected data
}
\description{
  Perform topographic correction using the Minnaert method.
  This code is modified from the code in the \code{landsat}
  package written by Sarah Goslee.  This version of the
  code has been altered from the \code{landsat} version to
  allow the option of using a sample of pixels for
  calculation of k in the Minnaert correction (useful when
  dealing with large images).
}
\details{
  See the help page for \code{minnaert} in the
  \code{landsat} package for additional details on the
  parameters.
}
\author{
  Sarah Goslee and Alex Zvoleff
}
\references{
  Sarah Goslee. Analyzing Remote Sensing Data in {R}: The
  {landsat} Package. Journal of Statistical Software, 2011,
  43:4, pg 1--25. http://www.jstatsoft.org/v43/i04/
}

