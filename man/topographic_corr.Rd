\name{topographic_corr}
\alias{topographic_corr}
\title{Topographically correct a raster}
\usage{
  topographic_corr(x, slopeaspect, sunelev, sunazimuth,
    method, filename = NULL, overwrite = FALSE,
    inparallel = FALSE, sampleindices = NULL)
}
\arguments{
  \item{x}{an image to correct}

  \item{slopeaspect}{a \code{RasterBrick} or
  \code{RasterStack} with two layers.  First layer should
  be the slope, second layer should be aspect. The slope
  and aspect are defined as in \code{slopeasp} in the
  \code{landsat} package.  \code{\link{slopeasp_par}} will
  output the slope and aspect using the proper definition
  and as a \code{RasterBrick}.}

  \item{sunelev}{sun elevation in degrees}

  \item{sunazimuth}{sun azimuth in degrees}

  \item{method}{the topographic correction method to use.
  See the help for \code{\link{topocorr}}.}

  \item{filename}{file on disk to save \code{Raster*} to
  (optional)}

  \item{overwrite}{whether to overwrite \code{filename} if
  it already exists}

  \item{inparallel}{whether to run correction in parallel
  using \code{foreach}}

  \item{sampleindices}{(optional) row-major indices of
  sample pixels to use in regression models used for some
  topographic correction methods (like Minnaert). Useful
  when handling very large images. See
  \code{\link{gridsample}} for one method of calculating
  these indices.}
}
\value{
  The topographically corrected image.
}
\description{
  Performs topographic correction using code based on
  \code{\link{topocorr}} from the \code{Landsat} package.
  The code in this package has been modifed from
  \code{\link{topocorr}} to allow using a subsample of the
  image for Minnaert k calculations, and to provide the
  option of running the topographic correction in parallel
  using \code{\link{foreach}}.
}
\examples{
# Mosaic the two ASTER DEM tiles needed to a Landsat image
DEM_mosaic <- mosaic(ASTER_V002_EAST, ASTER_V002_WEST, fun='mean')

# Crop and extend the DEM mosaic to match the Landsat image
matched_DEM <- match_rasters(L5TSR_1986, DEM_mosaic)
slopeaspect <- slopeasp_seq(matched_DEM)

# Apply the topographic correction
sunelev <- 90 - 44.97 # From metadata file
sunazimuth <- 124.37 # From metadata file
L5TSR_1986_topocorr <- topographic_corr(L5TSR_1986, slopeaspect, sunelev,
                                        sunazimuth, method='minslope')

plotRGB(L5TSR_1986, stretch='lin', r=3, g=2, b=1)
plotRGB(L5TSR_1986_topocorr, stretch='lin', r=3, g=2, b=1)
}

