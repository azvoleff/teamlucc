\name{extract_training_data}
\alias{extract_training_data}
\title{Extract training data for use in a classification}
\usage{
  extract_training_data(x, polys, classcol, testfrac = 0)
}
\arguments{
  \item{x}{a \code{Raster*} object from which training data
  will be extracted. Training data will be extracted from
  each layer in a \code{RasterBrick} or
  \code{RasterStack}.}

  \item{polys}{a \code{SpatialPolygonsDataFrame} with
  training data polygons}

  \item{classcol}{the name of the column containing the
  response variable (for example the land cover type of
  each pixel)}

  \item{testfrac}{a number from 0-1 giving the fraction of
  the data to be used for testing.}
}
\value{
  A data.frame with the training data. Each row will
  contain the response (the column chosen by
  \code{classcol}) as the first column, with the remaining
  columns containing the values at that location of each
  band in the raster stack.
}
\description{
  Extract training data for use in a classification
}
\examples{
L5TSR_1986 <- stack(system.file('extdata/L5TSR_1986.dat', package='teamr'))
data(L5TSR_1986_2001_training)
train_data <- extract_training_data(L5TSR_1986, L5TSR_1986_2001_training,
                                    "t1_class")
}

