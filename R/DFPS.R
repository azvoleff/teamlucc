#' Double-Window Flexible Pace Search (DFPS) threshold determination
#'
#' @export
#' @importFrom rgeos gDifference gBuffer
#' @param chg_polys \code{SpatialPolygonsDataFrame} with polygons of change 
#' areas surrounded by windows of no-change
#' @param chg_mag change magnitude \code{RasterLayer} from \code{CVAPS}
#' @param radius radius of no-change surrounding change area polygons
#' @param delta the minimum difference between Lmax and Lmin that will allow 
#' another pass through the search loop
#' @param m number of potential thresholds per search iteration (see Chen et 
#' al., 2003). Recommended to leave at default value.
#' @param maxiter maximum number of iterations of the main search process
#' @references Chen, J., P. Gong, C. He, R. Pu, and P. Shi. 2003.
#' Land-use/land-cover change detection using improved change-vector analysis.
#' Photogrammetric Engineering and Remote Sensing 69:369-380.
#' 
#' Chen, J., X. Chen, X. Cui, and J. Chen. 2011. Change vector analysis in
#' posterior probability space: a new method for land cover change detection.
#' IEEE Geoscience and Remote Sensing Letters 8:317-321.
DFPS <- function(chg_polys, chg_mag, radius=100, delta=.01, m=10, maxiter=20) {
    chg_pixels <- unlist(extract(chg_mag, chg_polys))
    nochg_polys <- gDifference(gBuffer(chg_polys, width=radius), chg_polys)
    nochg_pixels <- unlist(extract(chg_mag, nochg_polys))
    # Calculate total number of pixels (used later)
    A <- length(chg_pixels) + length(nochg_pixels)
    # Set initial values for Lmax and Lmin that ensure the below while loop 
    # will run.
    Lmax <- 1
    Lmin <- Lmax - 2 * delta
    min_threshold <- min(chg_pixels)
    max_threshold <- max(chg_pixels)
    n <- 0
    while ((Lmax - Lmin) > delta && n < maxiter) {
        p <- (max_threshold - min_threshold) / m
        thresholds <- seq(min_threshold, max_threshold, p)
        L <- c()
        for (threshold in thresholds) {
            A1 <- sum(chg_pixels > threshold)
            A2 <- sum(nochg_pixels > threshold)
            # Calculate A according to equation 4 in Chen et al. 2003
            L <- c(L, ((A1 - A2) * 100) / A)
        }
        kmax <- thresholds[match(max(L), L)]
        min_threshold <- kmax - p
        max_threshold <- kmax + p
        Lmin <- min(L)
        Lmax <- max(L)
        n <- n + 1
    }
    return(kmax)
}
