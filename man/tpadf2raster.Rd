\name{tpadf2raster}
\alias{tpadf2raster}
\title{Function to convert TIMESAT .tpa data.frame to an R raster.}
\usage{
  tpadf2raster(x, base_image, variable)
}
\arguments{
  \item{x}{A TPA data.frame as output by tpa2df}

  \item{base_image}{A string giving the location of a
  raster file to use for georeferencing the output raster.
  Use one of the original raster files that was input to
  TIMESAT.}

  \item{variable}{A string giving the variable name to
  write to a raster.  Can be one of: start, end, length,
  base_value, peak_time, peak_value, amp, left_deriv,
  right_deriv, large_integ, and small_integ.}
}
\value{
  A raster object
}
\description{
  Function to convert TIMESAT .tpa data.frame to an R
  raster.
}
\examples{
# TODO: Need to add examples here, and need to include a sample TIMESAT tpa
# file in the package data.
}

