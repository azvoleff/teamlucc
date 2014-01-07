\name{team_setup_dem}
\alias{team_setup_dem}
\title{Setup the DEM mosaic for a particular TEAM site}
\usage{
team_setup_dem(dem_list, sitecode, output_path, sample_image = NULL,
  n_cpus = 1)
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
}
\description{
Setup the DEM mosaic for a particular TEAM site
}
\examples{
\dontrun{
dem_list <- list('H:/Data/TEAM/VB/Rasters/DEM/ASTER/ASTGTM2_N09W084_dem.tif',
                 'H:/Data/TEAM/VB/Rasters/DEM/ASTER/ASTGTM2_N09W085_dem.tif',
                 'H:/Data/TEAM/VB/Rasters/DEM/ASTER/ASTGTM2_N10W084_dem.tif',
                 'H:/Data/TEAM/VB/Rasters/DEM/ASTER/ASTGTM2_N09W085_dem.tif')
team_setup_dem(dem_list, "VB", 'H:/Data/TEAM/VB/LCLUC_Analysis/')
}
}

