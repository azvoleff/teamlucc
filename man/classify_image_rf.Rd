\name{classify_image_rf}
\alias{classify_image_rf}
\title{Runs an image classification using a random forest classifier}
\usage{
classify_image_rf(x, train_data, class_probs = TRUE,
  use_training_flag = TRUE, train_control = NULL, tune_grid = NULL,
  use_rfe = FALSE, notify = print)
}
\arguments{
  \item{x}{a \code{Raster*} image with the predictor
  layer(s) for the classification}

  \item{train_data}{a \code{link{pixel_data}} object}

  \item{class_probs}{whether to also calculate and return
  the probabilities of membership for each class}

  \item{use_training_flag}{indicates whether to exclude
  data flagged as testing data when training the
  classifier. For this to work the input train_data
  \code{data.frame} must have a column named
  'training_flag' that indicates, for each pixel, whether
  that pixel is a training pixel (coded as TRUE) or testing
  pixel (coded as FALSE).}

  \item{train_control}{default is NULL (reasonable value
  will be set automatically).  For details see
  \code{\link{trainControl}}.}

  \item{tune_grid}{the training grid to be used for
  training the classifier. Must be a \code{data.frame} with
  two columns: ".sigma" and ".C".}

  \item{use_rfe}{whether to use Recursive Feature
  Extraction (RFE) as implemented in the \code{caret}
  package to select a subset of the input features to be
  used in the classification}

  \item{notify}{notifier to use (defaults to \code{print}
  function). See the \code{notifyR} package for one way of
  sending notifications from R. The \code{notify} function
  should accept a string as the only argument.}
}
\value{
a list with 3 elements: the trained classifier, the
predicted classes \code{RasterLayer} and the class
probabilities \code{RasterBrick}
}
\description{
Currently only supports classification using a random
forest
}
\details{
Processing can be done in parallel using all using the
cluster facilities in the \code{spatial.tools} package. To
enable clustering, call \code{beginCluster} before running
\code{classify_image_rf}.  To stop the cluster when
finished, call \code{endCluster}.
}
\examples{
\dontrun{
# Don't run long example
set.seed(0)
train_data_1986 <- get_pixels(L5TSR_1986, polys=L5TSR_1986_2001_training,
                              class_col="class_1986", training=.7)
classified_LT5SR_1986 <- classify_image_rf(L5TSR_1986, train_data_1986)

classified_LT5SR_1986$model
plot(classified_LT5SR_1986$pred_classes)
plot(classified_LT5SR_1986$pred_probs)
accuracy(classified_LT5SR_1986$model)
}
}

