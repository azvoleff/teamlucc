\name{topocorr_samp}
\alias{topocorr_samp}
\title{Topographic correction for satellite imagery}
\usage{
  topocorr_samp(x, slope, aspect, sunelev, sunazimuth,
    method = "cosine", na.value = NA, GRASS.aspect = FALSE,
    IL.epsilon = 1e-06, sampleindices = NULL)
}
\arguments{
  \item{x}{image as a \code{RasterLayer}}

  \item{slope}{the slope as a \code{RasterLayer}}

  \item{aspect}{the aspect as a \code{RasterLayer}}

  \item{sunelev}{sun elevation in degrees}

  \item{sunazimuth}{sun azimuth in degrees}

  \item{method}{the method to use for the topographic
  correction: cosine, improvedcosine, minnaert, minslope,
  ccorrection, gamma, SCS, or illumination}

  \item{na.value}{the value used to code no data values}

  \item{GRASS.aspect}{is aspect defined according to GRASS
  convetion}

  \item{IL.epsilon}{a small amount to add to calculated
  illumination values that are equal to zero to avoid
  division by zero resulting in Inf values}

  \item{sampleindices}{(optional) row-major indices of
  sample pixels to use in regression models used for some
  topographic correction methods (like Minnaert). Useful
  when handling very large images. See
  \code{\link{gridsample}} for one method of calculating
  these indices.}
}
\value{
  RasterBrick with two layers: 'slope' and 'aspect'
}
\description{
  Perform topographic correction using a number of
  different methods. This code is modified from the code in
  the \code{landsat} package by Sarah Goslee.  This version
  of the code has been altered from the \code{landsat}
  version to allow the option of using a sample of pixels
  for calculation of k in the Minnaert correction (useful
  when dealing with large images).
}
\details{
  See the help page for \code{topocorr} in the
  \code{landsat} package for details on the parameters.
}
\examples{
#TODO: add examples
}
\author{
  Sarah Goslee and Alex Zvoleff
}
\references{
  Sarah Goslee. Analyzing Remote Sensing Data in {R}: The
  {landsat} Package. Journal of Statistical Software, 2011,
  43:4, pg 1--25. http://www.jstatsoft.org/v43/i04/
}

