\name{accuracy}
\alias{accuracy}
\title{Calculate statistics on classification accuracy}
\usage{
  accuracy(model, test_data)
}
\arguments{
  \item{model}{a classification model with a \code{predict}
  method}

  \item{test_data}{a training/testing dataset as output by
  the \code{extract_training_data} function. If a
  'Training' column is included \code{accuracy} will use
  only the data not used in training for evaluating model
  accuracy.}
}
\value{
  list of accuracy statistics
}
\description{
  Used for calculating a contingency table with user's,
  producer's, and overall accuracies for an image
  classification.
}

