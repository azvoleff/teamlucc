\name{applyWindowed}
\alias{applyWindowed}
\title{Apply a raster function with edge effects over a series of blocks}
\usage{
applyWindowed(x, fun, out_file, edge = c(0, 0), ...)
}
\arguments{
  \item{x}{a \code{Raster*}}

  \item{fun}{the function to apply}

  \item{out_file}{path to save output raster to}

  \item{edge}{length 2 numberic with number of rows on top
  and bottom with edge effects, defined as c(top, bottom)}

  \item{...}{additional arguments to pass to \code{fun}}
}
\description{
This function can be useful when applying windowed
functions over a raster, as with \code{glcm}. This function
allows windows functions that have edge effects to be
applied over a raster in block-by-block fashion.
\code{applyWindowed} avoids the striping that would result
if the edge effects were ignored.
}
\examples{
\dontrun{
applyWindowed(L5TSR_1986, glcm, 'out.envi', edge=c(1, 3))
}
}

