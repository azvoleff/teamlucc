#' Topographically correct a raster
#'
#' Performs topographic correction using code based on \code{topocorr} from the 
#' \code{landsat} package by Sarah Goslee. The code in this package has been 
#' modifed from \code{topocorr} to allow using a subsample of the image for 
#' Minnaert k calculations, and to provide the option of running the 
#' topographic correction in parallel using \code{foreach}.
#'
#' This function will run in parallel if a parallel backend is registered with 
#' \code{\link{foreach}}.
#'
#' @export
#' @import foreach
#' @param x an image to correct
#' @param slopeaspect a \code{RasterBrick} or \code{RasterStack} with two 
#' layers.  The first layer should be the slope, the second layer should be 
#' the aspect. The slope and aspect are defined as in \code{terrain} in the 
#' \code{raster} package, and both should be in radians.
#' @param sunelev sun elevation in degrees
#' @param sunazimuth sun azimuth in degrees
#' @param method the topographic correction method to use. See the help for 
#' \code{topocorr} for more guidance on this.
#' @param sampleindices (optional) row-major indices of sample pixels to use in 
#' regression models used for some topographic correction methods (like 
#' Minnaert). Useful when handling very large images. See
#' \code{\link{gridsample}} for one method of calculating these indices.
#' @param scale_factor factor by which to multiply results. Useful if rounding 
#' results to integers (see \code{asinteger} argument).
#' @param asinteger whether to round results to nearest integer. Can be used to 
#' save space by saving results as, for example, an 'INT2S' \code{raster}.
#' @param DN_min minimum allowable pixel value after correction (values less 
#' than \code{DN_min} are set to NA)
#' @param DN_max maximum allowable pixel value after correction (values less 
#' than \code{DN_max} are set to NA)
#' @param ... additional arguments to pass to \code{minnaert_samp} or 
#' \code{topocorr_samp}, depending on chosen topographic correction method
#' @return The topographically corrected image as a \code{RasterLayer} or 
#' \code{RasterStack}
#' @references
#' Sarah Goslee. Analyzing Remote Sensing Data in {R}: The {landsat} Package.  
#' Journal of Statistical Software, 2011, 43:4, pg 1--25.  
#' http://www.jstatsoft.org/v43/i04/
#' @examples
#' \dontrun{
#' # Mosaic the two ASTER DEM tiles needed to a Landsat image
#' DEM_mosaic <- mosaic(ASTER_V002_EAST, ASTER_V002_WEST, fun='mean')
#' 
#' # Crop and extend the DEM mosaic to match the Landsat image
#' matched_DEM <- match_rasters(L5TSR_1986, DEM_mosaic)
#' slopeaspect <- terrain(matched_DEM, opt=c('slope', 'aspect')
#' 
#' # Apply the topographic correction
#' sunelev <- 90 - 44.97 # From metadata file
#' sunazimuth <- 124.37 # From metadata file
#' L5TSR_1986_topocorr <- topographic_corr(L5TSR_1986, slopeaspect, sunelev, 
#'                                         sunazimuth, method='minslope')
#' 
#' plotRGB(L5TSR_1986, stretch='lin', r=3, g=2, b=1)
#' plotRGB(L5TSR_1986_topocorr, stretch='lin', r=3, g=2, b=1)
#' }
topographic_corr <- function(x, slopeaspect, sunelev, sunazimuth, 
                             method='minnaert_full', sampleindices=NULL, 
                             scale_factor=1, asinteger=FALSE, DN_min=NULL, 
                             DN_max=NULL, ...) {
    if (!(class(x) %in% c('RasterLayer', 'RasterStack', 'RasterBrick'))) {
        stop('x must be a Raster* object')
    }
    if (!(class(slopeaspect) %in% c('RasterBrick', 'RasterStack'))) {
        stop('slopeaspect must be a RasterBrick or RasterStack object')
    }
    slope <- raster(slopeaspect, layer=1)
    aspect <- raster(slopeaspect, layer=2)
    stopifnot((sunelev >= 0) & (sunelev <= 90))
    stopifnot((sunazimuth >= 0) & (sunazimuth <= 360))

    # Set uncorr_layer to NULL to pass R CMD CHECK without notes
    uncorr_layer=NULL
    was_rasterlayer <- class(x) == "RasterLayer"
    # Convert x to stack so foreach will run
    x <- stack(x)
    corr_img <- foreach(uncorr_layer=unstack(x), .combine='addLayer', 
                        .multicombine=TRUE, .init=raster(), 
                        .packages=c('teamlucc', 'rgdal')) %dopar% {
        if (method == 'minnaert_full') {
            minnaert_data <- minnaert_samp(uncorr_layer, slope, aspect, 
                                           sunelev=sunelev, 
                                           sunazimuth=sunazimuth, 
                                           sampleindices=sampleindices, 
                                           DN_min=DN_min, DN_max=DN_max, ...)
            corr_layer <- minnaert_data$minnaert
        } else {
            corr_layer <- topocorr_samp(uncorr_layer, slope, aspect, 
                                        sunelev=sunelev, sunazimuth=sunazimuth, 
                                        method=method, 
                                        sampleindices=sampleindices,
                                        DN_min=DN_min, DN_max=DN_max, ...)
        }
    }
    if (was_rasterlayer) corr_img <- corr_img[[1]]
    names(corr_img) <- paste0(names(x), 'tc')
    if (scale_factor != 1) {
        corr_img <- corr_img * scale_factor
    }
    if (asinteger) {
        corr_img <- round(corr_img)
    }
    return(corr_img)
}
