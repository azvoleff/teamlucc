\name{gridsample}
\alias{gridsample}
\title{Draw a random sample from a grid laid out on a RasterLayer or matrix}
\usage{
  gridsample(x, horizcells = 10, vertcells = 10,
    nsamp = 10, rowmajor = FALSE, replace = FALSE)
}
\arguments{
  \item{x}{a matrix or RasterLayer to draw sample from}

  \item{horizcells}{how many cells to break the raster in
  horizontally (over the columns)}

  \item{vertcells}{how many cells to break the raster in
  vertically (over the rows)}

  \item{nsamp}{how many samples to draw from each grid
  cell}

  \item{rowmajor}{whether to return indices in row-major
  format (default is column-major). Row-major format is
  useful in conjunction with \code{Raster*} objects.}

  \item{replace}{whether to sample with replacement (within
  each grid cell)}
}
\value{
  vector of sample indices
}
\description{
  This function is used to subsample a \code{RasterLayer}
  or \code{matrix} by dividing the dataset into a grid of
  \code{horizcells} x \code{vertcells}, and by then drawing
  a sample of size \code{nsamp} from within each grid cell.
}
\note{
  TODO: Recode in C++ for speed.
}
\examples{
# Make a 100x100 matrix
x <- matrix(1:10000, nrow=100)
# Subsample the matrix by breaking it up into a 10 x 10 grid, and take 10
# random samples from each grid cell without replacement (these are the
# default parameters).
y <- gridsample(x)
}

