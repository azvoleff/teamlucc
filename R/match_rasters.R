#' Match the coordinate system and extent of two rasters
#'
#' @export
#' @param baseimg A /code{Raster*} to use as the base image. This layer will 
#' determine the output coordinate system.
#' @param matchimg A /code{Raster*} to match to the base image. If necessary 
#' the /code{matchimg} will be reprojected to match the coordinate system of 
#' the /code{baseimg}. The /code{matchimg} will then be cropped and extended to 
#' match the extent of the /code{baseimg}.
#' @param filename file on disk to save \code{Raster*} to (optional)
#' @param method the method to use if projection is needed to match image 
#' coordinate systems, or if resampling is needed to align image origins. Can 
#' be "ngb" for nearest-neighbor, or "binlinear" for bilinear interpolation
#' @param ... additional arguments to pass to \code{writeRaster} (such as 
#' datatype and filename)
#' @return The /code{matchimg} reprojected (if necessary), cropped, and 
#' extended to match the /code{baseimg}.
#' @details Note that \code{match_rasters} can run in parallel if 
#' \code{beginCluster()} is run prior to running \code{match_rasters}.
#' @examples
#' # Mosaic the two ASTER DEM tiles needed to a Landsat image
#' DEM_mosaic <- mosaic(ASTER_V002_EAST, ASTER_V002_WEST, fun='mean')
#' 
#' # Crop and extend the DEM mosaic to match the Landsat image
#' matched_DEM <- match_rasters(L5TSR_1986, DEM_mosaic)
match_rasters <- function(baseimg, matchimg, filename, method='bilinear',
                          ...) {
    if (projection(baseimg) != projection(matchimg)) {
        matchimg <- projectRaster(from=matchimg, to=baseimg, method=method)
    }

    # Crop overlapping area
    matchimg <- crop(matchimg, baseimg)
    if (any(res(matchimg) != res(baseimg))) {
        matchimg <- resample(matchimg, baseimg, method=method)
    }
    if (extent(matchimg) != extent(baseimg)) {
        matchimg <- extend(matchimg, baseimg)
    }

    if (!missing(filename)) {
            matchimg <- writeRaster(matchimg, filename=filename, ...)
    }
    return(matchimg)
}
