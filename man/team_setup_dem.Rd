\name{team_setup_dem}
\alias{team_setup_dem}
\title{Setup the DEM mosaic for a particular TEAM site}
\usage{
team_setup_dem(dem_path, sitecode, output_path, sample_image = NULL,
  n_cpus = 1, overwrite = FALSE, notify = print)
}
\arguments{
  \item{dem_list}{a list of digital elevation models (DEMs)
  that (when mosaiced) covers the full extent of all the
  images in the image_list.}

  \item{sitecode}{code to use as a prefix for all
  filenames}

  \item{output_path}{the path to use for the output}

  \item{sample_image}{a \code{Raster*} or the path to an
  image in a format readable by \code{raster}. The sample
  image will be used to set the projection system for the
  output DEM mosaic.}

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
Setup the DEM mosaic for a particular TEAM site
}
\examples{
\dontrun{
dem_path <- 'H:/Data/TEAM/VB/Rasters/DEM/ASTER'
team_setup_dem(dem_path, "VB", 'H:/Data/TEAM/VB/LCLUC_Analysis/')
}
}

