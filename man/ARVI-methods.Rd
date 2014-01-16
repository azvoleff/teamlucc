\docType{methods}
\name{ARVI}
\alias{ARVI}
\alias{ARVI,RasterLayer,RasterLayer,RasterLayer-method}
\alias{ARVI,matrix,matrix,matrix-method}
\alias{ARVI,numeric,numeric,numeric-method}
\title{Calculates the Atmospherically Resistant Vegetation Index (ARVI)}
\usage{
ARVI(blue, red, nir, ...)

\S4method{ARVI}{numeric,numeric,numeric}(blue, red, nir)

\S4method{ARVI}{matrix,matrix,matrix}(blue, red, nir)

\S4method{ARVI}{RasterLayer,RasterLayer,RasterLayer}(blue, red, nir, ...)
}
\arguments{
  \item{blue}{blue}

  \item{red}{red}

  \item{nir}{near-infrared}

  \item{...}{additional arguments as for
  \code{\link{writeRaster}}}
}
\description{
Calculates the Atmospherically Resistant Vegetation Index
(ARVI)
}
\examples{
ARVI_img <- ARVI(blue=raster(L5TSR_1986, layer=1), red=raster(L5TSR_1986,
                 layer=3), nir=raster(L5TSR_1986, layer=4))
plot(ARVI_img)
}
\references{
Kaufman, Y. J., and D. Tanre. 1996. Strategy for direct and
indirect methods for correcting the aerosol effect on
remote sensing: from AVHRR to EOS-MODIS. Remote Sensing of
Environment:65-79.
}

