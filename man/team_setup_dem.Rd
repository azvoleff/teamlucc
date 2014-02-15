\name{team_setup_dem}
\alias{team_setup_dem}
\title{Setup the DEM mosaic for a particular TEAM site}
\usage{
team_setup_dem(dem_path, output_path, aoi_file, n_cpus = 1,
  overwrite = FALSE, crop_to_aoi = FALSE, notify = print,
  verbose = FALSE)
}
\arguments{
  \item{dem_path}{a list of digital elevation models (DEMs)
  that (when mosaiced) covers the full extent of all the
  images in the image_list.}

  \item{aoi_file}{area of interest (AOI) shapefile to use
  as as bounding box when selecting DEMs, in a file format
  readable by \code{readOGR}}

  \item{output_path}{the path to use for the output}

  \item{n_cpus}{the number of CPUs to use for processes
  that can run in parallel}

  \item{overwrite}{whether to overwrite existing files
  (otherwise an error will be raised)}

  \item{crop_to_aoi}{whether to crop the dem to the
  supplied AOI, or to the Landsat path/row polygon for that
  particular path/row}

  \item{notify}{notifier to use (defaults to \code{print}
  function). See the \code{notifyR} package for one way of
  sending notifications from R. The \code{notify} function
  should accept a string as the only argument.}

  \item{verbose}{whether to print detailed status messages
  and timing information}
}
\description{
Setup the DEM mosaic for a particular TEAM site
}
\examples{
\dontrun{
dem_path <- 'H:/Data/TEAM/VB/Rasters/DEM/ASTER'
team_setup_dem(dem_path, "VB", 'H:/Data/TEAM/VB/LCLUC_Analysis/',
               list(c(15,53)))
}
}

