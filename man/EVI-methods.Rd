\docType{methods}
\name{EVI}
\alias{EVI}
\alias{EVI,RasterLayer,RasterLayer,RasterLayer-method}
\alias{EVI,matrix,matrix,matrix-method}
\alias{EVI,numeric,numeric,numeric-method}
\title{Calculates the Enhanced Vegetation Index (EVI)}
\usage{
EVI(blue, red, nir, ...)

\S4method{EVI}{numeric,numeric,numeric}(blue, red, nir)

\S4method{EVI}{matrix,matrix,matrix}(blue, red, nir)

\S4method{EVI}{RasterLayer,RasterLayer,RasterLayer}(blue, red, nir, ...)
}
\arguments{
  \item{blue}{blue}

  \item{red}{red}

  \item{nir}{near-infrared}

  \item{...}{additional arguments as for
  \code{\link{writeRaster}}}
}
\description{
Calculates the Enhanced Vegetation Index (EVI)
}
\examples{
EVI_img <- EVI(blue=raster(L5TSR_1986, layer=1), red=raster(L5TSR_1986, layer=3),
               nir=raster(L5TSR_1986, layer=4))
plot(EVI_img)
}
\references{
Huete, A. R., HuiQing Liu, and W. J. D. van Leeuwen. 1997.
The use of vegetation indices in forested regions: issues
of linearity and saturation. Pages 1966-1968 vol.4
Geoscience and Remote Sensing, 1997. IGARSS '97. Remote
Sensing - A Scientific Vision for Sustainable Development.,
1997 IEEE International.
}

