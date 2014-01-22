\docType{methods}
\name{scale_raster}
\alias{scale_raster}
\alias{scale_raster,RasterBrick,ANY-method}
\alias{scale_raster,RasterBrick-method}
\alias{scale_raster,RasterLayer,ANY-method}
\alias{scale_raster,RasterLayer-method}
\alias{scale_raster,RasterStack,ANY-method}
\alias{scale_raster,RasterStack-method}
\title{Scales a layer by power of a given integer and rounds to nearest integer}
\usage{
scale_raster(x, power_of = 10, max_out = 32767, round_output = TRUE,
  do_scaling = TRUE)

\S4method{scale_raster}{RasterLayer}(x, power_of, max_out, round_output,
  do_scaling)

\S4method{scale_raster}{RasterStack}(x, power_of, max_out, round_output,
  do_scaling)

\S4method{scale_raster}{RasterBrick}(x, power_of, max_out, round_output,
  do_scaling)
}
\arguments{
  \item{x}{a \code{Raster*} object}

  \item{power_of}{raster will be scaled using the highest
  possible power of this number}

  \item{max_out}{scale_factor will be chosen to ensure that
  the maximum and minimum values of the raster do not
  exceed this number}

  \item{do_scaling}{perform the scaling (if
  \code{do_scaling} is TRUE) or just return the scale
  factors (if \code{do_scaling} is FALSE)}

  \item{round_output}{whether to round the output to the
  nearest integer}
}
\description{
Useful for scaling a layer and rounding to integer so that
a layer can be saved as 'INT2' or 'INT2S'.
}

