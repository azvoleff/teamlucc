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
#' @param x time 0 posterior probability \code{Raster*}
#' @param y time 1 posterior probability \code{Raster*}
#' @param filename (optional) filename for output change magnitude 
#' \code{RasterLayer}
#' @param \code{Raster*} object with change direction image
#' @references Chen, J., P. Gong, C. He, R. Pu, and P. Shi. 2003.
#' Land-use/land-cover change detection using improved change-vector analysis.
#' Photogrammetric Engineering and Remote Sensing 69:369-380.
#' 
#' Chen, J., X. Chen, X. Cui, and J. Chen. 2011. Change vector analysis in 
#' posterior probability space: a new method for land cover change detection.  
#' IEEE Geoscience and Remote Sensing Letters 8:317-321.
chg_dir <- function(x, y, filename=NULL, ...) {
    if (proj4string(x) != proj4string(y)) {
        stop('Error: t0 and t1 coordinate systems do not match')
    }
    if (extent(x) != extent(y)) {
        stop('Error: t0 and t1 extents do not match')
    }
    if (nlayers(x) != nlayers(y)) {
        stop('Error: t0 and t1 probability maps have differing number of classes')
    }

    n_classes <- nlayers(x)
    if (n_classes < 2) {
        stop('Error: cannot calculate change probabilities for only one class')
    }

    # focal_hpc will only take one raster object as input, so stack the 
    # probability layers for time 1 and time 2 in a single RasterStack, then 
    # split them apart within the calc_chg_dir function:
    calc_chg_dir <- function(x, n_classes, ...) {
        t1p <- x[ , , seq(1, n_classes)]
        t2p <- x[ , , seq(n_classes + 1, dim(x)[3])]
        # Calculate change direction (eqns 5 and 6 in Chen 2011)
        dP <- array(t2p - t1p, dim=c(dim(x)[1], dim(x)[2], n_classes))
        #Eab <- dP %*% diag(n_classes)
        unit_vecs <- array(diag(n_classes), dim=c(3, 3))
        Eab <- apply(dP, c(1, 2), function(pixel) pixel %*% unit_vecs)
        chgdir <- apply(Eab, c(2, 3),
                        function(pixel) which(pixel == max(pixel)))
        chgdir <- array(chgdir, dim=c(dim(x)[1], dim(x)[2], 1))
        return(chgdir)
    }

    out <- focal_hpc(x=stack(x, y), fun=calc_chg_dir, 
                     args=list(n_classes=n_classes), filename=filename)

    return(out)
}
