\name{calc_texture}
\alias{calc_texture}
\title{Calculates a glcm texture for use in the glcm.R script}
\usage{
  calc_texture(Rast, statistics, base_indices,
    offset_indices, n_grey)
}
\arguments{
  \item{Rast}{a matrix containing the pixels to be used in
  the texture calculation}

  \item{statistics}{a list of strings naming the texture
  statistics to calculate}

  \item{base_indices}{a list of column-major indices giving
  the cells included in the base image}

  \item{offset_indices}{base_indices a list of column-major
  indices giving the cells ' included in the offset image}

  \item{n_grey}{number of grey levels to use in texture
  calculation}
}
\value{
  a list of length equal to the length of the
  \code{statistics} input parameter, containing the
  selected textures measures
}
\description{
  This function is called by the \code{\link{glcm}}
  function. It is not intended to be used directly.
}
\references{
  Sarah Goslee. Analyzing Remote Sensing Data in {R}: The
  {landsat} Package. Journal of Statistical Software, 2011,
  43:4, pg 1--25. http://www.jstatsoft.org/v43/i04/
}

