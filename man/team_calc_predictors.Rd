\name{team_calc_predictors}
\alias{team_calc_predictors}
\title{Calculate predictor layers for a classification}
\usage{
team_calc_predictors(image_dirs, output_path, n_cpus = 1, cleartmp = FALSE,
  overwrite = FALSE, notify = print)
}
\arguments{
  \item{image_dirs}{list of paths to a set of Landsat
  images that have been preprocessed by the
  \code{team_preprocess_landsat} function.}

  \item{output_path}{the path to use for the output}

  \item{n_cpus}{the number of CPUs to use for processes
  that can run in parallel}

  \item{cleartmp}{whether to clear temp files on each run
  through the loop}

  \item{overwrite}{whether to overwrite existing files
  (otherwise an error will be raised)}

  \item{notify}{notifier to use (defaults to \code{print}
  function). See the \code{notifyR} package for one way of
  sending notifications from R. The \code{notify} function
  should accept a string as the only argument.}
}
\description{
Calculate predictor layers for a classification
}

