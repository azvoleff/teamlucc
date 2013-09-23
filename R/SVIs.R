#' Calculates the Normalized Difference Vegetation Index (NDVI)
#'
#' @export
#' @param red \code{RasterLayer}
#' @param nir near-infrared \code{RasterLayer}
NDVI <- function(red, nir) {
  NDVI <- (nir - red) / (nir + red)
}

#' Calculates the Enhanced Vegetation Index (EVI)
#'
#' @export
#' @param blue \code{RasterLayer}
#' @param red \code{RasterLayer}
#' @param nir near-infrared \code{RasterLayer}
EVI <- function(blue, red, nir) {
  EVI <- 2.5*((nir - red) / (nir + 6.*red - 7.5*blue + 1.))
}

#' Calculates the Modified Soil-Adjusted Vegetation Index (MSAVI)
#'
#' @export
#' @param red \code{RasterLayer}
#' @param nir near-infrared \code{RasterLayer}
MSAVI2 <- function(red, nir) {
  MSAVI2 <- (2.*nir + 1 - sqrt((2.*nir + 1.)^2. - 8.*(nir - red)))/2.
}

#' Calculates the Atmospherically Resistant Vegetation Index (ARVI)
#'
#' @export
#' @param blue \code{RasterLayer}
#' @param red \code{RasterLayer}
#' @param nir near-infrared \code{RasterLayer}
ARVI <- function(blue, red, nir) {
    ARVI <- (nir - 2.*red - blue) / (nir + 2.*red - blue)
}
