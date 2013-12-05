\name{compcont}
\alias{compcont}
\title{Calculate a contingency table using the composite operator}
\usage{
  compcont(cont_t, ref_map = NULL, comp_map = NULL)
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
  contingency table
}
\description{
  This function calculates a cross tabulation for a map
  comparison using the composite operator recommended by
  Pontius and Cheuk (2006).
}
\references{
  Pontius, R. G., and M. Millones. 2011. Death to Kappa:
  birth of quantity disagreement and allocation
  disagreement for accuracy assessment. International
  Journal of Remote Sensing 32:4407-4429.
}

