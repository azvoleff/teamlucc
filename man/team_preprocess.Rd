\name{team_preprocess}
\alias{team_preprocess}
\title{Preprocess surface reflectance imagery from the Landsat CDR archive}
\usage{
team_preprocess(image_list, dem_list, n_cpus = 2, cleartmp = FALSE)
}
\arguments{
  \item{image_list}{a list of paths to Landsat CDR images
  that have been converted to bsq files with teampy.}

  \item{dem}{a list of digital elevation models (DEMs) that
  (when mosaiced) covers the full extent of all the images
  in the image_list.}

  \item{n_cpus}{the number of CPUs to use for processes
  that can run in parallel}

  \item{cleartmp}{whether to clear temp files on each run
  through the loop}
}
\description{
Preprocess surface reflectance imagery from the Landsat CDR
archive
}
\examples{
\dontrun{
image_list <- list('H:/Data/TEAM/VB/Rasters/Landsat/1986_037_LT5/proc/lndsr.LT50150531986037XXX18')
dem_list <- list('H:/Data/TEAM/VB/Rasters/DEM/ASTER/ASTGTM2_N09W084_dem.tif',
                 'H:/Data/TEAM/VB/Rasters/DEM/ASTER/ASTGTM2_N09W085_dem.tif',
                 'H:/Data/TEAM/VB/Rasters/DEM/ASTER/ASTGTM2_N10W084_dem.tif',
                 'H:/Data/TEAM/VB/Rasters/DEM/ASTER/ASTGTM2_N09W085_dem.tif')
team_preprocess(image_list, dem_list)
}
}

