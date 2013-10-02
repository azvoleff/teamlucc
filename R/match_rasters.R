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
#' @return The /code{matchimg} reprojected (if necessary), cropped, and 
#' extended to match the /code{baseimg}.
#' @details Note that if the \code{matchimg} needs to be reprojected,
#' \code{match_rasters} can run in parallel if \code{beginCluster()} is run 
#' prior to running \code{match_rasters}.
#' @examples
#' # Mosaic the four ASTER DEM tiles needed to cover the Landsat image
#' data(ASTER_V002_LL)
#' data(ASTER_V002_LR)
#' data(ASTER_V002_UR)
#' data(ASTER_V002_UL)
#' DEM_mosaic <- mosaic(ASTER_V002_LL, ASTER_V002_LR, ASTER_V002_UR, ASTER_V002_UL, fun='mean')
#' 
#' # Crop and extend the DEM mosaic to match the Landsat image
#' L5TSR_1986 <- raster(system.file('extdata/L5TSR_1986.dat', package='teamr'))
#' matched_DEM <- match_rasters(L5TSR_1986, DEM_mosaic)
match_rasters <- function(baseimg, matchimg, filename=NULL) {
    if (projection(baseimg) != projection(matchimg)) {
        message('Coordinate systems do not match - reprojecting matchimg...')
        matchimg <- projectRaster(matchimg, baseimg)
    }
    # First crop out any overlapping area
    message('Cropping matchimg to base...')
    outimg <- crop(matchimg, baseimg)
    # Now extend borders of cropped raster to match base raster
    message('Extending matchimg to base...')
    outimg <- extend(matchimg, outimg)
    #message('Resampling matchimg to base...')
    #resample(outimg, baseimg)
    if (!is.null(filename)) {
        writeRaster(outimg, filename)
    }
    return(outimg)
}
