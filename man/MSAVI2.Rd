\name{MSAVI2}
\alias{MSAVI2}
\title{Calculates the Modified Soil-Adjusted Vegetation Index (MSAVI)}
\usage{
  MSAVI2(red, nir)
}
\arguments{
  \item{red}{red \code{RasterLayer}}

  \item{nir}{near-infrared \code{RasterLayer}}
}
\description{
  Note that this avoids the need for calculating L by using
  the equation for MSAVI2 from Qi et al. (1994).
}
\examples{
MSAVI_img <- MSAVI2(red=raster(L5TSR_1986, layer=3),
                    nir=raster(L5TSR_1986, layer=4))
plot(MSAVI_img)
}
\references{
  Qi, J., A. Chehbouni, A. R. Huete, Y. H. Kerr, and S.
  Sorooshian. 1994. A modified soil adjusted vegetation
  index. Remote Sensing of Environment 48:119-126.
}

