#' Slope and aspect calculation for DEM \code{RasterLayer}
#'
#' This code calculate the slope and aspect for a DEM provided as a 
#' \code{RasterLayer}.
#'
#' The code for calculating the slope and aspect is taken directly from the 
#' \code{landsat} package by Sarah Goslee. See the help page for 
#' \code{slopeasp} in the \code{landsat} package for details on the parameters.
#'
#' @export
#' @import spatial.tools
#' @importFrom landsat slopeasp
#' @param DEM DEM as \code{RasterLayer}
#' @param EWkernel kernel to use for East-West gradient calculation
#' @param NSkernel kernel to use for North-South gradient calculation
#' @param smoothing positive integer for smoothing. 1 means no smoothing.
#' @param filename output file for 2-band slope and aspect layer stack 
#' (optional)
#' @param overwrite whether to overwrite output file if it exists
#' @return RasterBrick with two layers: 'slope' and 'aspect'
#' @author Sarah Goslee and Alex Zvoleff
#' @references
#' Sarah Goslee. Analyzing Remote Sensing Data in {R}: The {landsat} Package.  
#' Journal of Statistical Software, 2011, 43:4, pg 1--25.  
#' http://www.jstatsoft.org/v43/i04/
#' @examples
#' #TODO: Write examples
slopeasp_seq <- function(DEM, EWkernel, NSkernel, smoothing=1, filename='', 
                          overwrite=FALSE) {
    if (class(DEM) == 'SpatialGridDataFrame') {
        stop('DEM must be a Raster* object')
    }
    DEM_df <- as(DEM, "SpatialGridDataFrame")
    if (missing(EWkernel)) {
        EWkernel <- matrix(c(-1/8, 0, 1/8, -2/8, 0, 2/8, -1/8, 
            0, 1/8), ncol = 3, nrow = 3, byrow = TRUE)
    }
    if (missing(NSkernel)) {
        NSkernel <- matrix(c(1/8, 2/8, 1/8, 0, 0, 0, -1/8, -2/8, 
            -1/8), ncol = 3, nrow = 3, byrow = TRUE)
    }
    if (!identical(dim(NSkernel), dim(EWkernel))) {
        stop('NSkernel and EWkernel must have same dimensions')
    }
    slopeasp_spdf <- slopeasp(DEM_df, EWres=xres(DEM), EWkernel, 
                              NSres=yres(DEM), NSkernel=NSkernel, 
                              smoothing=smoothing)
    slope <- as(slopeasp_spdf$slope, 'RasterLayer')
    aspect <- as(slopeasp_spdf$aspect, 'RasterLayer')
    slopeasp_img <- brick(stack(slope, aspect), filename=filename, 
                          overwrite=overwrite)
    names(slopeasp_img) <- c('slope', 'aspect')
    return(slopeasp_img)
}
