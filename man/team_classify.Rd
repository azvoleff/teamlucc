\name{team_classify}
\alias{team_classify}
\title{Classify a preprocessed surface reflectance image}
\usage{
team_classify(predictor_file, train_shp, output_path, class_col = "Poly_Type",
  training = 0.6, n_cpus = 1, overwrite = FALSE, notify = print)
}
\arguments{
  \item{predictor_file}{a \code{Raster*} of predictor
  layers output by the \code{team_preprocess} function or
  path to an image stack in a format readable by the
  \code{raster} package.}

  \item{train_shp}{a training dataset as output by
  \code{extract_observed}}

  \item{output_path}{the path to use for the output}

  \item{class_col}{the name of the column containing the
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

