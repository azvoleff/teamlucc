\docType{methods}
\name{NDVI}
\alias{NDVI}
\alias{NDVI,RasterLayer,RasterLayer,RasterLayer-method}
\alias{NDVI,RasterLayer,RasterLayer-method}
\alias{NDVI,matrix,matrix,matrix-method}
\alias{NDVI,matrix,matrix-method}
\alias{NDVI,numeric,numeric,numeric-method}
\alias{NDVI,numeric,numeric-method}
\title{Calculates the Normalized Difference Vegetation Index (NDVI)}
\usage{
NDVI(red, nir, ...)

\S4method{NDVI}{numeric,numeric}(red, nir)

\S4method{NDVI}{matrix,matrix}(red, nir)

\S4method{NDVI}{RasterLayer,RasterLayer}(red, nir, ...)
}
\arguments{
  \item{red}{red}

  \item{nir}{near-infrared}

  \item{...}{additional arguments as for
  \code{\link{writeRaster}}}
}
\description{
Calculates the NDVI, defined as: (nir - red) / (nir + red).
}
\examples{
NDVI_img <- NDVI(red=raster(L5TSR_1986, layer=3), nir=raster(L5TSR_1986,
                layer=4))
plot(NDVI_img)
}

