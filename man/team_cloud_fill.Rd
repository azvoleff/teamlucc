\name{team_cloud_fill}
\alias{team_cloud_fill}
\title{Perform cloud fill for topographically corrected Landsat imagery}
\usage{
team_cloud_fill(data_dir, wrspath, wrsrow, start_date, end_date, output_path,
  base_date = NULL, n_cpus = 1, overwrite = FALSE, notify = print)
}
\arguments{
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
Uses the GNSPI algorithm from Xioalin Zhu. See
\code{\link{fill_clouds}} for details.
}

