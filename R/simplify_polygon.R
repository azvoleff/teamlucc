#' Count number of vertices in an sp polygon object
#'
#' @param poly_obj an sp polygon object
#' @return the number of vertices in the polygon
#' @examples
#' # TODO: Write an example.
nverts <- function(poly_obj) {
    n_verts <- sapply(poly_obj@polygons, function(y) 
                      nrow(y@Polygons[[1]]@coords))[[1]]
    if (is.null(n_verts)) {
        n_verts <- 0
    } else {
        # Need to subtract one as the starting coordinate and ending coordinate 
        # are identical, but appear twice - at beginning and at end of list.
        n_verts <- n_verts - 1
    }
    return(n_verts)
}

#' Simplify a polygon to contain less than a certain number of vertices
#'
#' Useful for simplifying area of interest (AOI) polygons for use with online 
#' data portals (like USGS EarthExplorer) that limit the number of vertices 
#' allowed in uploaded AOI shapefiles.
#'
#' @export
#' @importFrom rgeos gSimplify
#' @param poly_obj a polygon as an sp object
#' @param max_vertices the maximum number of vertices to allow in the 
#' simplified polygon
#' @param maxit the maximum number of iterations to use to simplify the polygon
#' @param multiplier a number used to increase the tolerance for each 
#' subsequent iteration, if the number of vertices in the simplified polygon is 
#' not less than /code{max_vertices}
#' @param initial_tolerance initial value for tolerance used to remove vertices 
#' from polygon.  If set to the default option, "dynamic", the code will 
#' automatically set the \code{initial_tolerance} to .01 * the length of the 
#' diagonal of the bounding box of the polygon. \code{initial_tolerance} can 
#' also be set to an arbitrary value, in the same units as the polygon object.
#' @return polygon with less than \code{max_vertices} vertices
#' @examples
#' # TODO: add an example
simplify_polygon <- function(poly_obj, max_vertices, maxit=100, 
                             multiplier=1.05, initial_tolerance='dynamic') {
    if (class(poly_obj) == 'SpatialPolygonsDataFrame') {
        poly_data <- poly_obj@data
    } else {
        poly_data <- NULL
    }
    n_parts <- sapply(poly_obj@polygons, function(x) length(x))
    if (length(n_parts) > 1)
        stop('poly_obj contains more than one polygon')
    else if (n_parts > 1)
        stop('poly_obj polygon is a multipart polygon')

    if (initial_tolerance == 'dynamic') {
        # Set the initial tolerance as the extent / 100
        ext <- extent(poly_obj)
        diag_seg_length <- sqrt((ext@xmax - ext@xmin)**2 +
                                (ext@ymax - ext@ymin)**2)
        tolerance <- diag_seg_length / 100
    } else {
        tolerance <- initial_tolerance
    }

    # Iterate, increasing tolerance, until polygon has less than maxpts 
    # vertices.
    n_verts <- nverts(poly_obj)
    n <- 0
    while ((n_verts > 0) && (n < maxit) && (n_verts > max_vertices)) {
        poly_obj <- gSimplify(poly_obj, tol=tolerance)
        n_verts <- nverts(poly_obj)
        tolerance <- tolerance * multiplier 
        n <- n + 1
    }
    if (n == maxit)
        warning(paste('Reached maximum iterations (', maxit, ')', sep=''))
    if (n_verts == 0)
        warning(paste('Simplified polygon has no vertices.'))
    else if (n_verts <= 2)
        warning(paste('Simplified polygon has only', n_verts, 'vertices.'))

    if (!is.null(poly_data)) {
        poly_obj <- SpatialPolygonsDataFrame(poly_obj, data=poly_data)
    }

    return(poly_obj)
}
