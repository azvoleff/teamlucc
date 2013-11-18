#' Calculates the Normalized Difference Vegetation Index (NDVI)
#'
#' @export
#' @param red red \code{RasterLayer}
#' @param nir near-infrared \code{RasterLayer}
#' @examples
#' L5TSR_1986 <- stack(system.file('extdata/L5TSR_1986.dat', package='teamr'))
#' NDVI_img <- NDVI(red=raster(L5TSR_1986, layer=3), nir=raster(L5TSR_1986, 
#'                 layer=4))
#' plot(NDVI_img)
NDVI <- function(red, nir) {
  NDVI <- (nir - red) / (nir + red)
}

#' Calculates the Enhanced Vegetation Index (EVI)
#'
#' @export
#' @param blue blue \code{RasterLayer}
#' @param red red \code{RasterLayer}
#' @param nir near-infrared \code{RasterLayer}
#' @references Huete, A. R., HuiQing Liu, and W. J. D. van Leeuwen. 1997. The
#' use of vegetation indices in forested regions: issues of linearity and
#' saturation. Pages 1966-1968 vol.4 Geoscience and Remote Sensing, 1997.
#' IGARSS '97. Remote Sensing - A Scientific Vision for Sustainable
#' Development., 1997 IEEE International.
#' @examples
#' L5TSR_1986 <- stack(system.file('extdata/L5TSR_1986.dat', package='teamr'))
#' EVI_img <- EVI(blue=raster(L5TSR_1986, layer=1), red=raster(L5TSR_1986, layer=3), 
#'                nir=raster(L5TSR_1986, layer=4))
#' plot(EVI_img)
EVI <- function(blue, red, nir) {
  EVI <- (2.5*(nir - red)) / (1 + nir + 6*red - 7.5*blue)
}

#' Calculates the Modified Soil-Adjusted Vegetation Index (MSAVI)
#'
#' @export
#' @param red red \code{RasterLayer}
#' @param nir near-infrared \code{RasterLayer}
#' @references Qi, J., A. Chehbouni, A. R. Huete, Y. H. Kerr, and S.
#' Sorooshian. 1994. A modified soil adjusted vegetation index. Remote Sensing
#' of Environment 48:119-126.
#' @examples
#' L5TSR_1986 <- stack(system.file('extdata/L5TSR_1986.dat', package='teamr'))
#' MSAVI2_img <- MSAVI2(red=raster(L5TSR_1986, layer=3), nir=raster(L5TSR_1986, 
#'                      layer=4))
#' plot(MSAVI2_img)
MSAVI2 <- function(red, nir) {
  MSAVI2 <- (2*nir + 1 - sqrt((2*nir + 1)^2 - 8*(nir - red)))/2
}

#' Calculates the Atmospherically Resistant Vegetation Index (ARVI)
#'
#' @export
#' @param blue blue \code{RasterLayer}
#' @param red red \code{RasterLayer}
#' @param nir near-infrared \code{RasterLayer}
#' @references Kaufman, Y. J., and D. Tanre. 1996. Strategy for direct and
#' indirect methods for correcting the aerosol effect on remote sensing: from
#' AVHRR to EOS-MODIS. Remote Sensing of Environment:65-79.
#' @examples
#' L5TSR_1986 <- stack(system.file('extdata/L5TSR_1986.dat', package='teamr'))
#' ARVI_img <- ARVI(blue=raster(L5TSR_1986, layer=1), red=raster(L5TSR_1986, 
#'                  layer=3), nir=raster(L5TSR_1986, layer=4))
#' plot(ARVI_img)
ARVI <- function(blue, red, nir) {
    ARVI <- (nir - 2*red - blue) / (nir + 2*red - blue)
}
