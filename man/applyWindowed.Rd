\name{applyWindowed}
\alias{applyWindowed}
\title{Apply a raster function with edge effects over a series of blocks}
\usage{
applyWindowed(x, fun, edge = c(0, 0), chunksize = NULL, ...)
}
\arguments{
  \item{x}{a \code{Raster*}}

  \item{fun}{the function to apply}

  \item{edge}{length 2 numberic with number of rows on top
  and bottom with edge effects, defined as c(top, bottom)}

  \item{chunksize}{the number of rows to read per block
  (passed to \code{raster} \code{blockSize} function.}

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
L5TSR_1986_b1 <- raster(L5TSR_1986, layer=1)
min_x <- cellStats(L5TSR_1986_b1, 'min')
max_x <- cellStats(L5TSR_1986_b1, 'max')
applyWindowed(L5TSR_1986_b1, glcm, edge=c(1, 3), min_x=min_x, max_x=max_x)
}
}

