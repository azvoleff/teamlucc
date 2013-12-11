\name{classify_image}
\alias{classify_image}
\title{Runs an image classification}
\usage{
  classify_image(x, train_data,
    pred_classes_filename = NULL,
    pred_probs_filename = NULL, train_grid = NULL,
    classProbs = FALSE, use_training_flag = TRUE)
}
\arguments{
  \item{x}{a \code{Raster*} image with the predictor
  layer(s) for the classification}

  \item{train_data}{a data table with a column labeled 'y'
  with the observed classes, and one or more columns with
  the values of predictor(s) at each location.}

  \item{pred_classes_filename}{the filename to use to save
  the predicted classes \code{RasterLayer}. If 'NULL' (the
  default) a temporary file will be used if necessary (if
  the size of the output raster exceeds available memory).}

  \item{pred_probs_filename}{the filename to use to save
  the class probabilities \code{RasterLayer} or
  \code{RasterBrick}. If 'NULL' (the default) a temporary
  file will be used if necessary (if the size of the output
  raster exceeds available memory).}

  \item{train_grid}{the training grid to be used for
  training the classifier. Must be a \code{data.frame} with
  two columns: ".sigma" and ".C".}

  \item{classProbs}{whether to also calculate and return
  the probabilities of membership for each class}

  \item{use_training_flag}{indicates whether to exclude
  data flagged as testing data when training the
  classifier. For this to work the input train_data
  \code{data.frame} must have a column named 'Training'
  that indicates, for each pixel, whether that pixel is a
  training pixel (coded as TRUE) or testing pixel (coded as
  FALSE).}
}
\value{
  a list with 3 elements: the trained classifier, the
  predicted classes \code{RasterLayer} and the class
  probabilities \code{RasterBrick}
}
\description{
  Currently only supports classification using a support
  vector machine (SVM).
}
\details{
  processing can be done in parallel using all available
  CPUs through the use of the cluster facilities in the
  \code{spatial.tools} package. To enable clustering, call
  \code{sfQuickInit} before running \code{classify_image}.
  To stop the cluster when finished, call
  \code{sfQuickStop}.
}
\examples{
\dontrun{
# Don't run long example
# Load example datasets
L5TSR_1986 <- stack(system.file('extdata/L5TSR_1986.dat', package='teamr'))
train_data_1986 <- extract_training_data(L5TSR_1986,
                                         polys=L5TSR_1986_2001_training,
                                         classcol="t1_class", training=.7)
# Supply a small training grid for classify_image to save processing time for
# the purposes of this example - in normal use, train_grid can be left
# unspecified.
classified_LT5SR_1986 <- classify_image(L5TSR_1986, train_data_1986,
    train_grid=data.frame(.sigma=.0495, .C=0.5), classProbs=TRUE)

classified_LT5SR_1986$model
plot(classified_LT5SR_1986$pred_classes)
plot(classified_LT5SR_1986$pred_probs)
}
}

