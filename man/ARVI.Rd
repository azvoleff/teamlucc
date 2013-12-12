\name{ARVI}
\alias{ARVI}
\title{Calculates the Atmospherically Resistant Vegetation Index (ARVI)}
\usage{
  ARVI(blue, red, nir)
}
\arguments{
  \item{blue}{blue \code{RasterLayer}}

  \item{red}{red \code{RasterLayer}}

  \item{nir}{near-infrared \code{RasterLayer}}
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
  Kaufman, Y. J., and D. Tanre. 1996. Strategy for direct
  and indirect methods for correcting the aerosol effect on
  remote sensing: from AVHRR to EOS-MODIS. Remote Sensing
  of Environment:65-79.
}

