#' Generate a SpatialPolygonDataFrame of raster extents
#'
#' Also includes the filename associated with each raster object. Useful for 
#' providing the \code{dem_extents} argument to the 
#' \code{\link{auto_setup_dem}} function.
#'
#' @export
#' @importFrom maptools spRbind
#' @param rast_list a \code{Raster*} object, or \code{list} of\code{Raster*} objects
#' @return \code{SpatialPolygonDataFrame} with the extent of each raster object 
#' as a polygon, with a "filename" attribute giving the filename for the raster 
#' object from with each extent is derived.
get_extent_polys <- function(rast_list) {
    if (!is.list(rast_list)) rast_list <- list(rast_list)

    proj4strings <- lapply(rast_list, function(x) proj4string(x))
    if (!all(proj4strings == proj4strings[[1]])) {
        stop('every raster in rast_list must have the same projection')
    }

    extents <- lapply(rast_list, function(x) extent(x))
    filenames <- lapply(rast_list, function(x) filename(x))
    # Convert extents to a list of SpatialPolygons objects
    extent_sps_list <- lapply(extents, function(x) as(x, 'SpatialPolygons'))

    # Convert from list of SpatialPolygons objects to a single SpatialPolygons 
    # object
    extent_sps <- extent_sps_list[[1]]
    if (length(extent_sps_list) > 1) {
        for (n in 2:length(extent_sps_list)) {
            extent_sps <- spRbind(extent_sps, spChFIDs(extent_sps_list[[n]], 
                                                 as.character(n)))
        }
    }

    # Finally convert the SpatialPolygons object into a 
    # SpatialPolygonsDataFrame that also includes the filename of the raster 
    # associated with each extent polygon as an attribute
    extent_polys <- SpatialPolygonsDataFrame(extent_sps, 
                                             data=data.frame(filename=unlist(filenames)))
    proj4string(extent_polys) <- proj4strings[[1]]

    extent_polys$filename <- as.character(extent_polys$filename)

    return(extent_polys)
}
