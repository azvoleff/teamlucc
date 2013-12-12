\name{training_sample}
\alias{training_sample}
\title{Generate random sample polygons from a raster layer}
\usage{
  training_sample(x, size, side = 2 * xres(x),
    fields = c(), out_file = NULL, validate = TRUE,
    validate.layer = 1, exp = 2)
}
\arguments{
  \item{x}{a \code{Raster*}}

  \item{size}{the sample size (number of sample polygons to
  return)}

  \item{side}{the length each side of the sample polygon}

  \item{fields}{a list of fields to include in the output
  \code{SpatialPolygonsDataFrame}}

  \item{out_file}{(optional) shapefile to save to save
  output}

  \item{validate}{whether to check that all sample polygons
  lie within image area (defined as having no NAs within
  the polygon)}

  \item{validate.layer}{the index within \code{x} of the
  layer that should be used for validation}

  \item{exp}{multiplier used to draw larger initial sample
  to account for the loss of sample polygons during
  validation. Increase this value if the final sample have
  fewer sample polygons than desired \code{size}}
}
\value{
  a \code{SpatialPolygonsDataFrame} with \code{n} polygons
}
\description{
  Useful for gathering training data for an image
  classification. With the default settings, the output
  polygons will be perfectly aligned with the pixels in the
  input raster.
}
\details{
  The \code{validate} parameter allows excluding no data
  areas from the output polygons (though the simple method
  this function uses to enforce this means that pixels
  within the border areas of the image will have a slightly
  lower probability of selection than those in the
  interior).
}
\examples{
\dontrun{
set.seed(0)
L5TSR_1986_b1 <- raster(system.file('extdata/L5TSR_1986.dat',
                        package='teamr'))
training_polys <- training_sample(L5TSR_1986_b1, 10,
                                  side=6*xres(L5TSR_1986_b1))
plot(L5TSR_1986_b1)
plot(training_polys, add=TRUE)
}
}

