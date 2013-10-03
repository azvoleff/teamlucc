#' Slope and aspect calculation in parallel for DEM \code{RasterLayer}
#'
#' This code calculate the slope and aspect for a DEM provided as a 
#' \code{RasterLayer}. The computations will be done in parallel if 
#' \code{sfQuickInit()} is called before \code{slopeasp_par}.
#'
#' The code for calculating the slope and aspect is taken directly from the 
#' \code{landsat} package by Sarah Goslee. See the help page for 
#' \code{slopeasp} in the \code{landsat} package for details on the parameters.
#'
#' @export
#' @param x DEM as \code{RasterLayer}
#' @param EWkernel kernel to use for East-West gradient calculation
#' @param NSkernel kernel to use for North-South gradient calculation
#' @param smoothing positive integer for smoothing. 1 means no smoothing.
#' @param filename output file for 2-band slope and aspect layer stack 
#' (optional)
#' @return list with two elements: 'slope' (a \code{RasteLayer} with the 
#' calculated slope) and 'aspect',  (a \code{RasteLayer} with the calculated 
#' aspect).
#' @references
#' Sarah Goslee. Analyzing Remote Sensing Data in {R}: The {landsat} Package.  
#' Journal of Statistical Software, 2011, 43:4, pg 1--25.  
#' http://www.jstatsoft.org/v43/i04/
#' @examples
#' #TODO: Write examples
slopeasp_par <- function(x, EWkernel, NSkernel, smoothing=1, filename=NULL) {
    EWres <- xres(x)
    NSres <- yres(x)
    if (missing(EWkernel)) {
        EWkernel <- matrix(c(-1/8, 0, 1/8, -2/8, 0, 2/8, -1/8, 
            0, 1/8), ncol = 3, nrow = 3, byrow = TRUE)
    }
    if (missing(NSkernel)) {
        NSkernel <- matrix(c(1/8, 2/8, 1/8, 0, 0, 0, -1/8, -2/8, 
            -1/8), ncol = 3, nrow = 3, byrow = TRUE)
    }
    # Convert EWkernel and NSkernel to 3d arrays for use with focal_hpc.
    EWkernel <- array(EWkernel, dim=c(dim(EWkernel)[1], dim(EWkernel)[2], 1))
    NSkernel <- array(NSkernel, dim=c(dim(NSkernel)[1], dim(NSkernel)[2], 1))
    if (!identical(dim(NSkernel), dim(EWkernel))) {
        stop('NSkernal and EWkernel must have same dimensions')
    }
    # Make function to pass to focal_hpc for calculating EW.mat and NW.mat.  
    # Return a 2 band layer stack: EW.mat as band 1, NS.mat as band 2.
    apply_kernel <- function(x, EWkernel, EWres, NSkernel, NSres, ...) {
        EW.mat <- sum(x * EWkernel) / EWres
        NS.mat <- sum(x * NSkernel) / NSres
        return(c(EW.mat, NS.mat))
    }
    print('got here 0')
    grad.mat <- focal_hpc(x, fun=apply_kernel,
                          args=list(EWkernel, EWres, NSkernel, NSres),
                          window_dims=dim(NSkernel))
    # Make function to pass to focal_hpc for calculating slope and aspect.  
    # Result will be returned as a layer stack, with slope in band 1, aspect in 
    # band 2.
    calc_slope <- function(x, smoothing, ...) {
        EW.mat <- x[ , , 1]
        NS.mat <- x[ , , 2]
        slope <- atan(sqrt(EW.mat^2 + NS.mat^2)/smoothing)
        slope <- (180/pi) * slope
        aspect <- 180 - (180/pi) * atan(NS.mat/EW.mat) + 90 * (EW.mat/abs(EW.mat))
        aspect[slope == 0] <- 0
        return(array(c(slope, aspect), dim=c(dim(x)[1], dim(x)[2], 2)))
    }
    print('got here 1')
    slopeasp_img <- focal_hpc(x=grad.mat, fun=calc_slope,
                        args=list(smoothing=smoothing), filename=filename)
    print('got here 2')
    return(list(slope=raster(slopeasp_img, layer=1),
           aspect=raster(slopeasp_img, layer=2)))
}

