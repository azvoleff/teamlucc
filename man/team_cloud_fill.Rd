\name{team_cloud_fill}
\alias{team_cloud_fill}
\title{Perform cloud fill for Landsat imagery}
\usage{
team_cloud_fill(data_dir, wrspath, wrsrow, start_date, end_date,
  base_date = NULL, fast = FALSE, n_cpus = 1, notify = print,
  threshold = 1, max_iter = 5, ...)
}
\arguments{
  \item{data_dir}{folder where input images are located,
  with filenames as output by the
  \code{\link{team_preprocess_landsat}} function}

  \item{wrspath}{World Reference System (WRS) path}

  \item{wrsrow}{World Reference System (WRS) row}

  \item{start_date}{start date of period from which images
  will be chosen to fill cloudy areas in the base image (as
  \code{Date} object)}

  \item{end_date}{end date of period from which images will
  be chosen to fill cloudy areas in the the base image (as
  \code{Date} object)}

  \item{base_date}{ideal date for base image (base image
  will be chosen as the image among the available images
  that is closest to this date). If NULL, then the base
  image will be the image with the lowest cloud cover.}

  \item{fast}{if \code{TRUE}, use the CLOUD_REMOVE_FAST.pro
  script. If \code{FALSE}, use the CLOUD_REMOVE.pro
  script.}

  \item{n_cpus}{the number of CPUs to use for processes
  that can run in parallel}

  \item{threshold}{maximum percent cloud cover allowable in
  base image (cloud fill iterate until percent cloud cover
  in base image is below this value) or until more than
  \code{max_iter} iterations have been run}

  \item{max_iter}{maximum number of times to run cloud fill
  script}

  \item{notify}{notifier to use (defaults to \code{print}
  function).  See the \code{notifyR} package for one way of
  sending notifications from R.  The \code{notify} function
  should accept a string as the only argument.}

  \item{...}{additional arguments passed to
  \code{\link{cloud_remove}}, such as \code{DN_min},
  \code{DN_max}, \code{use_IDL}, \code{fast},
  \code{verbose}, etc. See \code{\link{cloud_remove}}.}
}
\value{
\code{Raster*} object with cloud filled image.
}
\description{
Uses an R/C++ implementation of the NSPI algorithm from
Xioalin Zhu. See \code{\link{cloud_remove}} for details. In
hilly areas, cloud fill should be done after topographic
correction.
}

