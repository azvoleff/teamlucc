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
This function preprocesses surface reflectance imagery from
the Landsat Climate Data Record (CDR) archive. This
function can reproject CDR tiles to match the projection of
a given \code{aoi}, crop the tiles to match the \code{aoi}
or a common WRS-2 path/row polygon, mask missing data and
clouds out of the CDR tiles, and perform topographic
correction.
}
\details{
Prior to running \code{team_preprocess_landsat}, the
\code{\link{pespa_extract}} and should be used to extract
the original zipfiles supplied by USGS.  The
\code{\link{unstack_ledapscdr}} function should then be
used to unstack the HDF format files and convert them into
ENVI binary format. To perform topographic correction with
\code{team_preprocess_landsat}, first run the
\code{\link{team_setup_dem}} function to preprocess a set
of DEM tiles. Then run \code{team_preprocess_landsat} with
the \code{tc=TRUE} option.
}
\seealso{
\code{\link{espa_extract}}, \code{\link{unstack_ledapscdr},
\code{\link{team_setup_dem}}
}

