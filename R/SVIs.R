#' Calculates the Normalized Difference Vegetation Index (NDVI)
#'
#' Calculates the NDVI, defined as: (nir - red) / (nir + red).
#'
#' @export
#' @docType methods
#' @rdname NDVI-methods
#' @param red red
#' @param nir near-infrared
#' @param ... additional arguments as for \code{\link{writeRaster}}
#' @examples
#' NDVI_img <- NDVI(red=raster(L5TSR_1986, layer=3), nir=raster(L5TSR_1986, 
#'                 layer=4))
#' plot(NDVI_img)
setGeneric("NDVI", function(red, nir, ...) {
    standardGeneric("NDVI")
})

NDVI_calc <- function(red, nir) {
    v <- (nir - red) / (nir + red)
    return(v)
}

#' @rdname NDVI-methods
#' @aliases NDVI,numeric,numeric,numeric-method
setMethod("NDVI", signature(red="numeric", nir="numeric"),
    function(red, nir) {
        NDVI_calc(red, nir)
    }
)

#' @rdname NDVI-methods
#' @aliases NDVI,matrix,matrix,matrix-method
setMethod("NDVI", signature(red="matrix", nir="matrix"),
    function(red, nir) {
        ret <- NDVI_calc(red, nir)
        return(ret)
    }
)

#' @rdname NDVI-methods
#' @aliases NDVI,RasterLayer,RasterLayer,RasterLayer-method
#' @importFrom raster overlay
setMethod("NDVI", signature(red="RasterLayer", nir="RasterLayer"),
    function(red, nir, ...) {
        ret  <- overlay(red, nir, fun=function(red, nir) {
                NDVI_calc(red, nir)
            } , ...)
        return(ret)
    }
)

#' Calculates the Enhanced Vegetation Index (EVI)
#'
#' @export
#' @docType methods
#' @rdname EVI-methods
#' @importFrom raster overlay
#' @param blue blue
#' @param red red
#' @param nir near-infrared
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
setGeneric("EVI", function(blue, red, nir, ...) {
    standardGeneric("EVI")
})

EVI_calc <- function(blue, red, nir) {
    v <- (2.5*(nir - red)) / (1 + nir + 6*red - 7.5*blue)
    return(v)
}

#' @rdname EVI-methods
#' @aliases EVI,numeric,numeric,numeric-method
setMethod("EVI", signature(blue="numeric", red="numeric", nir="numeric"),
    function(blue, red, nir) {
        EVI_calc(blue, red, nir)
    }
)

#' @rdname EVI-methods
#' @aliases EVI,numeric,numeric,numeric-method
setMethod("EVI", signature(blue="matrix", red="matrix", nir="matrix"),
    function(blue, red, nir) {
        ret <- EVI_calc(blue, red, nir)
        return(ret)
    }
)

#' @rdname EVI-methods
#' @aliases EVI,RasterLayer,RasterLayer,RasterLayer-method
#' @importFrom raster overlay
setMethod("EVI", signature(blue="RasterLayer", red="RasterLayer", 
                           nir="RasterLayer"),
    function(blue, red, nir, ...) {
        ret <- overlay(blue, red, nir, fun=function(blue, red, nir) {
                EVI_calc(blue, red, nir)
            }, ...)
        return(ret)
    }
)

#' Calculates the Modified Soil-Adjusted Vegetation Index (MSAVI)
#'
#' Note that this avoids the need for calculating L by using the equation for 
#' MSAVI2 from Qi et al. (1994).
#'
#' @export
#' @docType methods
#' @rdname MSAVI2-methods
#' @importFrom raster overlay
#' @param red red
#' @param nir near-infrared
#' @param ... additional arguments as for \code{\link{writeRaster}}
#' @references Qi, J., A. Chehbouni, A. R. Huete, Y. H. Kerr, and S.
#' Sorooshian. 1994. A modified soil adjusted vegetation index. Remote Sensing
#' of Environment 48:119-126.
#' @examples
#' MSAVI_img <- MSAVI2(red=raster(L5TSR_1986, layer=3),
#'                     nir=raster(L5TSR_1986, layer=4))
#' plot(MSAVI_img)
setGeneric("MSAVI2", function(red, nir, ...) {
    standardGeneric("MSAVI2")
})

MSAVI2_calc <- function(red, nir) {
    v <- (2*nir + 1 - sqrt((2*nir + 1)^2 - 8*(nir - red)))/2
    return(v)
}

#' @rdname MSAVI2-methods
#' @aliases MSAVI2,numeric,numeric,numeric-method
setMethod("MSAVI2", signature(red="numeric", nir="numeric"),
    function(red, nir) {
        MSAVI2_calc(red, nir)
    }
)

#' @rdname MSAVI2-methods
#' @aliases MSAVI2,matrix,matrix,matrix-method
setMethod("MSAVI2", signature(red="matrix", nir="matrix"),
    function(red, nir) {
        ret <- MSAVI2_calc(red, nir)
        return(ret)
    }
)

#' @rdname MSAVI2-methods
#' @aliases MSAVI2,RasterLayer,RasterLayer,RasterLayer-method
#' @importFrom raster overlay
setMethod("MSAVI2", signature(red="RasterLayer", nir="RasterLayer"),
    function(red, nir, ...) {
        ret  <- overlay(red, nir, fun=function(red, nir) {
                MSAVI2_calc(red, nir)
            } , ...)
        return(ret)
    }
)

#' Calculates the Atmospherically Resistant Vegetation Index (ARVI)
#'
#' @export
#' @docType methods
#' @rdname ARVI-methods
#' @importFrom raster overlay
#' @param blue blue
#' @param red red
#' @param nir near-infrared
#' @param ... additional arguments as for \code{\link{writeRaster}}
#' @references Kaufman, Y. J., and D. Tanre. 1996. Strategy for direct and
#' indirect methods for correcting the aerosol effect on remote sensing: from
#' AVHRR to EOS-MODIS. Remote Sensing of Environment:65-79.
#' @examples
#' ARVI_img <- ARVI(blue=raster(L5TSR_1986, layer=1), red=raster(L5TSR_1986, 
#'                  layer=3), nir=raster(L5TSR_1986, layer=4))
#' plot(ARVI_img)
setGeneric("ARVI", function(blue, red, nir, ...) {
    standardGeneric("ARVI")
})


ARVI_calc <- function(blue, red, nir) {
    v <- (nir - 2*red - blue) / (nir + 2*red - blue)
    return(v)
}

#' @rdname ARVI-methods
#' @aliases ARVI,numeric,numeric,numeric-method
setMethod("ARVI", signature(blue="numeric", red="numeric", nir="numeric"),
    function(blue, red, nir) {
        ARVI_calc(blue, red, nir)
    }
)

#' @rdname ARVI-methods
#' @aliases ARVI,matrix,matrix,matrix-method
setMethod("ARVI", signature(blue="matrix", red="matrix", nir="matrix"),
    function(blue, red, nir) {
        ret <- ARVI_calc(blue, red, nir)
        return(ret)
    }
)

#' @rdname ARVI-methods
#' @aliases ARVI,RasterLayer,RasterLayer,RasterLayer-method
#' @importFrom raster overlay
setMethod("ARVI", signature(blue="RasterLayer", red="RasterLayer", 
                            nir="RasterLayer"),
    function(blue, red, nir, ...) {
        ret <- overlay(blue, red, nir, fun=function(blue, red, nir) {
                ARVI_calc(blue, red, nir)
            }, ...)
        return(ret)
    }
)
