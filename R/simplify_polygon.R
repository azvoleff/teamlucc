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

#' Simplify a TEAM site ZOI polygon for use with EarthExplorer
#'
#' @export
#' @param ZOI_poly a ZOI polygon as an sp object
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
#' @return The TEAM site ZOI boundary polygon simplified for use with online 
#' data portals that limit the number of vertices allowed in uploaded AOI 
#' shapefiles.
#' @examples
#' # TODO: Write an example. Will need to add a polygon object to the package.
#' # ZOI_poly <- readOGR('H:/Data/TEAM/VB/Vectors', 'VB_ZOI_GEO')
#' # simplify_ZOI(ZOI_poly, 3)
simplify_ZOI <- function(ZOI_poly, max_vertices, maxit=100, multiplier=1.25, 
                         initial_tolerance='dynamic') {
    n_parts <- sapply(ZOI_poly@polygons, function(x) length(x))
    if (length(n_parts) > 1)
        stop('Error: ZOI shapefile contains more than one polygon')
    else if (n_parts > 1)
        stop('Error: ZOI polygon is a multipart polygon')

    if (initial_tolerance == 'dynamic') {
        # Set the initial tolerance as the extent / 100
        ext <- extent(ZOI_poly)
        diag_seg_length <- sqrt((ext@xmax - ext@xmin)**2 +
                                (ext@ymax - ext@ymin)**2)
        tolerance <- diag_seg_length / 100
    } else if (!is.numeric(tolerance) | (length(tolerance) != 1)) {
        stop('Error: tolerance must be "dynamic" or a numeric of length one')
    } else {
        tolerance <- initial_tolerance
    }

    # Iterate, increasing tolerance, until polygon has less than maxpts 
    # vertices.
    n_verts <- nverts(ZOI_poly)
    n <- 0
    while ((n_verts > 0) && (n < maxit) && (n_verts > max_vertices)) {
        ZOI_poly <- gSimplify(ZOI_poly, tol=tolerance)
        n_verts <- nverts(ZOI_poly)
        tolerance <- tolerance * multiplier 
        n <- n + 1
    }
    if (n == maxit)
        warning(paste('Reached maximum iterations (', maxit, ')', sep=''))
    if (n_verts == 0)
        warning(paste('Simplified polygon has no vertices.'))
    else if (n_verts <= 2)
        warning(paste('Simplified polygon has only', n_verts, 'vertices.'))

    return (ZOI_poly)
}

