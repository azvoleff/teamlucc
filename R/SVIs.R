#' Calculates the Normalized Difference Vegetation Index (NDVI)
#'
#' @export
#' @param red red \code{RasterLayer}
#' @param nir near-infrared \code{RasterLayer}
#' @param asinteger whether to round results to nearest integer
#' @param ... additional arguments as for \code{\link{writeRaster}}
#' @examples
#' NDVI_img <- NDVI(red=raster(L5TSR_1986, layer=3), nir=raster(L5TSR_1986, 
#'                 layer=4))
#' plot(NDVI_img)
NDVI <- function(red, nir, asinteger=FALSE, ...) {
    ret  <- overlay(red, nir, function(red, nir) {
        res <- (nir - red) / (nir + red)
        if(asinteger) {
            res <- round(res)
        }
        return(res)
    }, ...)
    names(ret) <- 'NDVI'
    return(ret)
}

#' Calculates the Enhanced Vegetation Index (EVI)
#'
#' @export
#' @param blue blue \code{RasterLayer}
#' @param red red \code{RasterLayer}
#' @param nir near-infrared \code{RasterLayer}
#' @param asinteger whether to round results to nearest integer
#' @param ... additional arguments as for \code{\link{writeRaster}}
#' @references Huete, A. R., HuiQing Liu, and W. J. D. van Leeuwen. 1997. The
#' use of vegetation indices in forested regions: issues of linearity and
#' saturation. Pages 1966-1968 vol.4 Geoscience and Remote Sensing, 1997.
#' IGARSS '97. Remote Sensing - A Scientific Vision for Sustainable
#' Development., 1997 IEEE International.
#' @examples
#' EVI_img <- EVI(blue=raster(L5TSR_1986, layer=1), red=raster(L5TSR_1986, layer=3), 
#'                nir=raster(L5TSR_1986, layer=4))
#' plot(EVI_img)
EVI <- function(blue, red, nir, asinteger=FALSE, ...) {
    ret <- overlay(blue, red, nir, function(blue, red, nir) {
        res <- (2.5*(nir - red)) / (1 + nir + 6*red - 7.5*blue)
        if(asinteger) {
            res <- round(res)
        }
        return(res)
    }, ...)
    names(ret) <- 'EVI'
    return(ret)
}

#' Calculates the Modified Soil-Adjusted Vegetation Index (MSAVI)
#'
#' Note that this avoids the need for calculating L by using the equation for 
#' MSAVI2 from Qi et al. (1994).
#'
#' @export
#' @param red red \code{RasterLayer}
#' @param nir near-infrared \code{RasterLayer}
#' @param asinteger whether to round results to nearest integer
#' @param ... additional arguments as for \code{\link{writeRaster}}
#' @references Qi, J., A. Chehbouni, A. R. Huete, Y. H. Kerr, and S.
#' Sorooshian. 1994. A modified soil adjusted vegetation index. Remote Sensing
#' of Environment 48:119-126.
#' @examples
#' MSAVI_img <- MSAVI2(red=raster(L5TSR_1986, layer=3),
#'                     nir=raster(L5TSR_1986, layer=4))
#' plot(MSAVI_img)
MSAVI2 <- function(red, nir, asinteger=FALSE, ...) {
    ret <- overlay(red, nir, function(red, nir) {
        res <- 2*nir + 1 - sqrt((2*nir + 1)^2 - 8*(nir - red)))/2
        if(asinteger) {
            res <- round(res)
        }
        return(res)
    }, ...)
    names(ret) <- 'MSAVI2'
    return(ret)
}

#' Calculates the Atmospherically Resistant Vegetation Index (ARVI)
#'
#' @export
#' @param blue blue \code{RasterLayer}
#' @param red red \code{RasterLayer}
#' @param nir near-infrared \code{RasterLayer}
#' @param asinteger whether to round results to nearest integer
#' @param ... additional arguments as for \code{\link{writeRaster}}
#' @references Kaufman, Y. J., and D. Tanre. 1996. Strategy for direct and
#' indirect methods for correcting the aerosol effect on remote sensing: from
#' AVHRR to EOS-MODIS. Remote Sensing of Environment:65-79.
#' @examples
#' ARVI_img <- ARVI(blue=raster(L5TSR_1986, layer=1), red=raster(L5TSR_1986, 
#'                  layer=3), nir=raster(L5TSR_1986, layer=4))
#' plot(ARVI_img)
ARVI <- function(blue, red, nir, asinteger=FALSE, ...) {
    ret <- overlay(blue, red, nir, function(blue, red, nir) {
        res <- (nir - 2*red - blue) / (nir + 2*red - blue)
        if(asinteger) {
            res <- round(res)
        }
        return(res)
    }, ...)
    names(ret) <- 'ARVI'
    return(ret)
}
