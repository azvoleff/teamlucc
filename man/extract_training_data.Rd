\name{extract_training_data}
\alias{extract_training_data}
\title{Extract training data for use in a classification}
\usage{
  extract_training_data(x, polys, classcol, training = 1)
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

  \item{training}{indicator of which polygons to use in
  training. Can be: 1) a string giving the name of a column
  indicating whether each polygon is to be used in training
  (column equal to TRUE) or in testing (column equal to
  FALSE), or 2) a logical vector of length equal to
  length(polys), or 3) a number between 0 and 1 indicating
  the fraction of the polygons to be randomly selected for
  use in training.}
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
train_data <- extract_training_data(L5TSR_1986, L5TSR_1986_2001_training,
"class_1986")
}

