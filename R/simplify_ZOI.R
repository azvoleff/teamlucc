#' Simplify a TEAM site ZOI polygon for use with EarthExplorer
#'
#' @export
#' @param shapefile a shapefile of a TEAM site ZOI boundary polygon
#' @return The TEAM site ZOI boundary polygon simplified for use with online 
#' data portals that limit the number of vertices allowed in uploaded AOI 
#' shapefiles.
#' @examples
#' # TODO: Write an example.
#' require(rgdal)
#' require(rgeos)
#' require(sp)
ZOI_poly <- readOGR('H:/Data/TEAM/VB/Vectors', 'VB_ZOI_GEO')
simplify_ZOI(ZOI_poly, 10, 3)
# TODO: Might need to load 'shapefiles' package to get extent function

simplify_ZOI <- function(ZOI_poly, max_vertices, maxit=100) {
    n_parts <- sapply(ZOI_poly@polygons, function(x) length(x))
    if (length(n_parts) > 1)
        stop('Error: ZOI shapefile contains more than one polygon')
    if (n_parts > 1)
        stop('Error: ZOI polygon is a multipart polygon')

    # Set the initial tolerance as the extent / 100
    ext <- extent(ZOI_poly)
    diag_seg_length <- sqrt((ext@xmax - ext@xmin)**2 +
                            (ext@ymax - ext@ymin)**2)
    tolerance <- diag_seg_length / 100

    # Iterate, increasing tolerance, until polygon has less than maxpts 
    # vertices.
    n_verts <- sapply(ZOI_poly@polygons, function(y) nrow(y@Polygons[[1]]@coords))
    n <- 0
    while ((n_verts > max_vertices) & (n < maxit)) {
        ZOI_poly <- gSimplify(ZOI_poly, tol=tolerance)
        n_verts <- sapply(ZOI_poly@polygons, function(y) nrow(y@Polygons[[1]]@coords))
        tolerance <- tolerance * 1.5
        print(tolerance)
        print(n_verts)
        n <- n + 1
    }
    if (n == maxit)
        warning(paste('reached maximum iterations (', maxit, ')', sep=''))
    return (ZOI_poly)
}

