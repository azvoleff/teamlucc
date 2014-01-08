\name{team_classify}
\alias{team_classify}
\title{Classify a preprocessed surface reflectance image}
\usage{
team_classify(predictors, train_data, output_path, n_cpus, overwrite = FALSE,
  notify = print)
}
\arguments{
  \item{predictors}{a \code{Raster*} of predictor layers
  output by the \code{team_preprocess} function or path to
  an image stack in a format readable by the \code{raster}
  package.}

  \item{train_data}{a training dataset as output by
  \code{extract_training_data}}

  \item{output_path}{the path to use for the output}

  \item{n_cpus}{the number of CPUs to use for processes
  that can run in parallel}

  \item{overwrite}{whether to overwrite existing files
  (otherwise an error will be raised)}

  \item{notify}{notifier to use (defaults to \code{print}
  function). See the \code{notifyR} package for one way of
  sending notifications from R. The \code{notify} function
  should accept a string as the only argument.}
}
\description{
First the image should be preprocesed using the
\code{team_preprocess} function.
}
\examples{
#TODO: Add example
}

