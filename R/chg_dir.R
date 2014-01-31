#' Change Direction Image for CVAPS
#'
#' This code calculate the change direction image for the Change Vector 
#' Analysis in Posterior Probability Space (CVAPS) method of Chen et al. 2011.  
#' Use the change direction image in conjunction with the change magnitude 
#' image from \code{chg_dir}, and \code{DFPS} to use the Double Window Flexible 
#' Pace Search method (Chen et al. 2003) to determine the threshold to use to 
#' map areas of change and no-change.
#'
#' @export
#' @importFrom spatial.tools rasterEngine
#' @param t1p time 0 posterior probability \code{Raster*}
#' @param t2p time 1 posterior probability \code{Raster*}
#' @param filename (optional) filename for output change direction
#' \code{RasterLayer}
#' @param ... additional parameters to pass to rasterEngine
#' @return \code{Raster*} object with change direction image
#' @details Processing can be done in parallel using all using the cluster 
#' facilities in the \code{spatial.tools} package. To enable clustering, call 
#' \code{beginCluster} before running \code{classify_image}.  To stop the 
#' cluster when finished, call \code{endCluster}.
#' @references Chen, J., P. Gong, C. He, R. Pu, and P. Shi. 2003.
#' Land-use/land-cover change detection using improved change-vector analysis.
#' Photogrammetric Engineering and Remote Sensing 69:369-380.
#' 
#' Chen, J., X. Chen, X. Cui, and J. Chen. 2011. Change vector analysis in 
#' posterior probability space: a new method for land cover change detection.  
#' IEEE Geoscience and Remote Sensing Letters 8:317-321.
chg_dir <- function(t1p, t2p, filename=NULL, ...) {
    if (proj4string(t1p) != proj4string(t2p)) {
        stop('t0 and t1 coordinate systems do not match')
    }
    if (extent(t1p) != extent(t2p)) {
        stop('t0 and t1 extents do not match')
    }
    if (nlayers(t1p) != nlayers(t2p)) {
        stop('t0 and t1 probability maps have differing number of classes')
    }

    n_classes <- nlayers(t1p)
    if (n_classes == 1) {
        stop('cannot calculate change probabilities for only one class')
    }

    calc_chg_dir <- function(t1p, t2p, n_classes, ...) {
        # Calculate change direction (eqns 5 and 6 in Chen 2011)
        dP <- array(t2p - t1p, dim=c(dim(t1p)[1], dim(t1p)[2], n_classes))
        unit_vecs <- array(diag(n_classes), dim=c(n_classes, n_classes))
        Eab <- apply(dP, c(1, 2), function(pixel) pixel %*% unit_vecs)
        chgdir <- apply(Eab, c(2, 3),
                        function(pixel) which(pixel == max(pixel)))
        chgdir <- array(chgdir, dim=c(dim(t1p)[1], dim(t1p)[2], 1))
        return(chgdir)
    }
    out <- rasterEngine(t1p=t1p, t2p=t2p, fun=calc_chg_dir, 
                        args=list(n_classes=n_classes), filename=filename,
                        outbands=1, ...)

    return(out)
}
