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
#' @param overwrite whether to overwrite \code{filename} if it already exists
#' @return The /code{matchimg} reprojected (if necessary), cropped, and 
#' extended to match the /code{baseimg}.
#' @details Note that if the \code{matchimg} needs to be reprojected,
#' \code{match_rasters} can run in parallel if \code{beginCluster()} is run 
#' prior to running \code{match_rasters}.
#' @examples
#' # Mosaic the two ASTER DEM tiles needed to a Landsat image
#' DEM_mosaic <- mosaic(ASTER_V002_EAST, ASTER_V002_WEST, fun='mean')
#' 
#' # Crop and extend the DEM mosaic to match the Landsat image
#' matched_DEM <- match_rasters(L5TSR_1986, DEM_mosaic)
match_rasters <- function(baseimg, matchimg, filename='', overwrite=FALSE) {
    if (projection(baseimg) != projection(matchimg)) {
        matchimg <- projectRaster(matchimg, baseimg)
    }
    # First crop out any overlapping area
    outimg <- crop(matchimg, baseimg)
    # Now extend borders of cropped raster to match base raster
    if (extent(outimg) != extent(baseimg)) {
        outimg <- extend(outimg, baseimg, filename=filename, 
                         overwrite=overwrite, datatype=dataType(outimg))
    } else if (filename != '') {
        outimg <- writeRaster(outimg, filename=filename, overwrite=overwrite, 
                              datatype=dataType(outimg))
    }
    return(outimg)
}
