\name{mixing_model}
\alias{mixing_model}
\title{Fractional mixing model implemented with support vector machine (SVM)}
\usage{
  mixing_model(x, train_data, filename = NULL)
}
\arguments{
  \item{x}{a \code{Raster*} image with the predictor
  layer(s) for the classification}

  \item{train_data}{a data table with a column labeled 'y'
  with the observed classes, and one or more columns with
  the values of predictor(s) at each location.}

  \item{filename}{(optional) filename for output image}
}
\value{
  a \code{RasterBrick} with the proportion of each pixel
  within each of the classes listed in the training data
}
\description{
  TODO: Write docs.
}
\references{
  Brown, M., H. G.  Lewis, and S.  R.  Gunn.  2000.  Linear
  spectral mixture models and support vector machines for
  remote sensing. Geoscience and Remote Sensing, IEEE
  Transactions on 38:2346-2360.
}

