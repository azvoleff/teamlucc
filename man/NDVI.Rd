\name{NDVI}
\alias{NDVI}
\title{Calculates the Normalized Difference Vegetation Index (NDVI)}
\usage{
  NDVI(red, nir)
}
\arguments{
  \item{red}{red \code{RasterLayer}}

  \item{nir}{near-infrared \code{RasterLayer}}
}
\description{
  Calculates the Normalized Difference Vegetation Index
  (NDVI)
}
\examples{
L5TSR_1986 <- stack(system.file('extdata/L5TSR_1986.dat', package='teamr'))
NDVI_img <- NDVI(red=raster(L5TSR_1986, layer=3), nir=raster(L5TSR_1986,
                layer=4))
plot(NDVI_img)
}

