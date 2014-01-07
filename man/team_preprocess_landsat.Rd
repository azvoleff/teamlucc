\name{team_preprocess_landsat}
\alias{team_preprocess_landsat}
\title{Preprocess surface reflectance imagery from the Landsat CDR archive}
\usage{
team_preprocess_landsat(image_list, dem, slopeaspect, sitecode, output_path,
  n_cpus = 1, cleartmp = FALSE, overwrite = FALSE)
}
\arguments{
  \item{image_list}{a list of paths to Landsat CDR images
  that have been converted to bsq files with teampy.}

  \item{dem}{path to a digital elevation model (DEM)
  covering the full extent of all the images in
  \code{image_list}. See \code{team_setup_dem} for a
  function simplifying this.}

  \item{slopeaspect}{path to a two layer raster stack with
  slope and aspect calculated from the above DEM}

  \item{sitecode}{code to use as a prefix for all
  filenames}

  \item{output_path}{the path to use for the output}

  \item{n_cpus}{the number of CPUs to use for processes
  that can run in parallel}

  \item{cleartmp}{whether to clear temp files on each run
  through the loop}

  \item{overwrite}{whether to overwrite existing files
  (otherwise an error will be raised)}
}
\description{
Preprocess surface reflectance imagery from the Landsat CDR
archive
}
\examples{
\dontrun{
image_list <- c('H:/Data/TEAM/VB/Rasters/Landsat/1986_037_LT5/proc/lndsr.LT50150531986037XXX18.bsq',
                'H:/Data/TEAM/VB/Rasters/Landsat/2001_014_LT5/proc/lndsr.LT50150532001014AAA01.bsq',
                'H:/Data/TEAM/VB/Rasters/Landsat/2012_021_LE7/proc/lndsr.LE70150532012021EDC00.bsq')
dem <- 'H:/Data/TEAM/VB/LCLUC_Analysis/VB_dem_mosaic.envi'
slopeaspect <- 'H:/Data/TEAM/VB/LCLUC_Analysis/VB_dem_mosaic_slopeaspect.envi'
team_preprocess(image_list, dem, slopeaspect, "VB",
'H:/Data/TEAM/VB/LCLUC_Analysis', 3, TRUE)
}
}

