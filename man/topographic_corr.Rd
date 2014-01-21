\name{topographic_corr}
\alias{topographic_corr}
\title{Topographically correct a raster}
\usage{
topographic_corr(x, slopeaspect, sunelev, sunazimuth, method, filename = "",
  inparallel = FALSE, sampleindices = NULL, scale_factor = 1,
  asinteger = FALSE, ...)
}
\arguments{
  \item{x}{an image to correct}

  \item{slopeaspect}{a \code{RasterBrick} or
  \code{RasterStack} with two layers.  The first layer
  should be the slope, the second layer should be the
  aspect. The slope and aspect are defined as in
  \code{terrain} in the \code{raster} package.}

  \item{sunelev}{sun elevation in degrees}

  \item{sunazimuth}{sun azimuth in degrees}

  \item{method}{the topographic correction method to use.
  See the help for \code{\link{topocorr}} for more guidance
  on this.}

  \item{filename}{file on disk to save \code{Raster*} to
  (optional)}

  \item{inparallel}{whether to run correction in parallel
  using \code{foreach}}

  \item{sampleindices}{(optional) row-major indices of
  sample pixels to use in regression models used for some
  topographic correction methods (like Minnaert). Useful
  when handling very large images. See
  \code{\link{gridsample}} for one method of calculating
  these indices.}

  \item{scale_factor}{factor by which to multiply results.
  Useful if rounding results to integers (see
  \code{asinteger} argument).}

  \item{asinteger}{whether to round results to nearest
  integer. Can be used to save space by saving results as,
  for example, an 'INT2S' \code{raster}.}

  \item{...}{additional arguments to pass to
  \code{writeRaster} (such as datatype and filename)}
}
\value{
The topographically corrected image as a \code{RasterLayer}
or \code{RasterStack}
}
\description{
Performs topographic correction using code based on
\code{topocorr} from the \code{Landsat} package. The code
in this package has been modifed from \code{topocorr} to
allow using a subsample of the image for Minnaert k
calculations, and to provide the option of running the
topographic correction in parallel using
\code{\link{foreach}}.
}
\examples{
\dontrun{
# Mosaic the two ASTER DEM tiles needed to a Landsat image
DEM_mosaic <- mosaic(ASTER_V002_EAST, ASTER_V002_WEST, fun='mean')

# Crop and extend the DEM mosaic to match the Landsat image
matched_DEM <- match_rasters(L5TSR_1986, DEM_mosaic)
slopeaspect <- terrain(matched_DEM, opt=c('slope', 'aspect')

# Apply the topographic correction
sunelev <- 90 - 44.97 # From metadata file
sunazimuth <- 124.37 # From metadata file
L5TSR_1986_topocorr <- topographic_corr(L5TSR_1986, slopeaspect, sunelev,
                                        sunazimuth, method='minslope')

plotRGB(L5TSR_1986, stretch='lin', r=3, g=2, b=1)
plotRGB(L5TSR_1986_topocorr, stretch='lin', r=3, g=2, b=1)
}
}

