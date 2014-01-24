\name{classify_image}
\alias{classify_image}
\title{Runs an image classification}
\usage{
classify_image(x, train_data, pred_classes_filename = NULL,
  pred_probs_filename = NULL, classProbs = TRUE, use_training_flag = TRUE,
  tune_length = 8, train_control = NULL, tune_grid = NULL)
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

  \item{classProbs}{whether to also calculate and return
  the probabilities of membership for each class}

  \item{use_training_flag}{indicates whether to exclude
  data flagged as testing data when training the
  classifier. For this to work the input train_data
  \code{data.frame} must have a column named 'Training'
  that indicates, for each pixel, whether that pixel is a
  training pixel (coded as TRUE) or testing pixel (coded as
  FALSE).}

  \item{tune_length}{the number of levels of each parameter
  that should be tried by \code{train} when training the
  classifier}

  \item{train_control}{default is NULL (reasonable value
  will be set automatically).  For details see
  \code{\link{trainControl}}.}

  \item{train_grid}{the training grid to be used for
  training the classifier. Must be a \code{data.frame} with
  two columns: ".sigma" and ".C".}
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
processing can be done in parallel using all available CPUs
through the use of the cluster facilities in the
\code{spatial.tools} package. To enable clustering, call
\code{sfQuickInit} before running \code{classify_image}. To
stop the cluster when finished, call \code{sfQuickStop}.
}
\examples{
\dontrun{
# Don't run long example
set.seed(0)
train_data_1986 <- extract_training_data(L5TSR_1986,
                                         polys=L5TSR_1986_2001_training,
                                         classcol="class_1986", training=.7)
classified_LT5SR_1986 <- classify_image(L5TSR_1986, train_data_1986)

classified_LT5SR_1986$model
plot(classified_LT5SR_1986$pred_classes)
plot(classified_LT5SR_1986$pred_probs)
summary(accuracy(classified_LT5SR_1986$model))
}
}

