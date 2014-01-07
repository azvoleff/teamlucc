\name{team_classify}
\alias{team_classify}
\title{Classify a preprocessed surface reflectance image}
\usage{
team_classify(predictors, train_data)
}
\arguments{
  \item{predictors}{a \code{Raster*} of predictor layers
  output by the \code{team_preprocess} function or path to
  an image stack in a format readable by the \code{raster}
  package.}

  \item{train_data}{a training dataset as output by
  \code{extract_training_data}}
}
\description{
First the image should be preprocesed using the
\code{team_preprocess} function.
}
\examples{
#TODO: Add example
}

