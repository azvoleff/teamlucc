\name{team_preprocess_landsat}
\alias{team_preprocess_landsat}
\title{Preprocess surface reflectance imagery from the Landsat CDR archive}
\usage{
team_preprocess_landsat(image_dirs, prefix, tc = FALSE, aoi = NULL,
  dem_path = NULL, output_path = NULL, n_cpus = 1, cleartmp = FALSE,
  overwrite = FALSE, notify = print, verbose = FALSE)
}
\arguments{
  \item{image_dirs}{list of paths to a set of Landsat CDR
  image files in ENVI format as output by the
  \code{unstack_ledapscdr} function.}

  \item{output_path}{the path to use for the output
  (optional - if NULL then output images will be saved
  alongside the input images in the same folder).}

  \item{prefix}{string to use as a prefix for all
  filenames}

  \item{tc}{whether to topographically correct imagery}

  \item{aoi}{area of interest (AOI), as a
  \code{SpatialPolygonsDataFrame}, to use as as bounding
  box when selecting DEMs}

  \item{dem_path}{path to a set of DEMs as output by
  \code{team_setup_dem} (only required if tc=TRUE)}

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

  \item{verbose}{whether to print detailed status messages
  and timing information}
}
\description{
Preprocess surface reflectance imagery from the Landsat CDR
archive
}
\examples{
\dontrun{
image_dirs <- c('H:/Data/TEAM/VB/Rasters/Landsat/1986_037_LT5/proc',
                'H:/Data/TEAM/VB/Rasters/Landsat/2001_014_LT5/proc',
                'H:/Data/TEAM/VB/Rasters/Landsat/2012_021_LE7/proc')
team_preprocess(image_dirs, 'H:/Data/TEAM/VB/LCLUC_Analysis', 'VB'
'H:/Data/TEAM/VB/LCLUC_Analysis', n_cpus=3)
}
}

