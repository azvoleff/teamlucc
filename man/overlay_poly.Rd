\name{overlay_poly}
\alias{overlay_poly}
\title{Make plot of image with overlaid polygon}
\usage{
overlay_poly(x, y, out, width = 6.5, height = 6.5, dpi = 300, ...)
}
\arguments{
  \item{x}{image as a \code{Raster*} object}

  \item{out}{filename for output image. The extension of
  this file will determine the output file format (png,
  pdf, etc.).}

  \item{width}{width (in inches) of output image}

  \item{height}{height (in inches) of output image}

  \item{dpi}{DPI for output image}

  \item{...}{additional arguments to \code{ggsave}}
}
\value{
used for writing plot to a file
}
\description{
Useful for quick plots showing overlap between an area of
interest (AOI) and a satellite image.
}

