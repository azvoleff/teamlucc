\name{match_rasters}
\alias{match_rasters}
\title{Match the coordinate system and extent of two rasters}
\usage{
  match_rasters(baseimg, matchimg, filename = NULL,
    overwrite = FALSE)
}
\arguments{
  \item{baseimg}{A /code{Raster*} to use as the base image.
  This layer will determine the output coordinate system.}

  \item{matchimg}{A /code{Raster*} to match to the base
  image. If necessary the /code{matchimg} will be
  reprojected to match the coordinate system of the
  /code{baseimg}. The /code{matchimg} will then be cropped
  and extended to match the extent of the /code{baseimg}.}

  \item{filename}{file on disk to save \code{Raster*} to
  (optional)}

  \item{overwrite}{whether to overwrite \code{filename} if
  it already exists}
}
\value{
  The /code{matchimg} reprojected (if necessary), cropped,
  and extended to match the /code{baseimg}.
}
\description{
  Match the coordinate system and extent of two rasters
}
\details{
  Note that if the \code{matchimg} needs to be reprojected,
  \code{match_rasters} can run in parallel if
  \code{beginCluster()} is run prior to running
  \code{match_rasters}.
}
\examples{
# Mosaic the two ASTER DEM tiles needed to a Landsat image
data(ASTER_V002_EAST, ASTER_V002_WEST)
DEM_mosaic <- mosaic(ASTER_V002_EAST, ASTER_V002_WEST, fun='mean')

# Crop and extend the DEM mosaic to match the Landsat image
L5TSR_1986 <- stack(system.file('extdata/L5TSR_1986.dat', package='teamr'))
matched_DEM <- match_rasters(L5TSR_1986, DEM_mosaic)
}

