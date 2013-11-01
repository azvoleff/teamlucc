#' Topographically correct a raster
#'
#' Performs topographic correction using \code{\link{topocorr}} from the 
#' \code{Landsat} package.
#'
#' @export
#' @import foreach
#' @importFrom landsat topocorr minnaert
#' @param x An image to correct.
#' @param sunelev Sun elevation in degrees.
#' @param sunazimuth Sun azimuth in degrees.
#' @param slopeaspect \code{RasterBrick} or \code{RasterStack} with two layers. 
#' First layer should be the slope, second layer should be aspect. The slope 
#' and aspect are defined as in \code{slopeasp} in the \code{landsat} package.  
#' \code{\link{slopeasp_par}} will output the slope and aspect using the proper 
#' definition and as a \code{RasterBrick}.
#' @param method the topographic correction method to use. See the help for 
#' \code{\link{topocorr}}.
#' @param filename file on disk to save \code{Raster*} to (optional)
#' @param overwrite whether to overwrite \code{filename} if it already exists
#' @param inparallel whether to run correction in parallel using \code{foreach}
#' @param usesample whether to subsample the data with \code{gridsample} 
#' (useful when handling very large images)
#' @param ... Additional arguments to be passed to \code{topocorr_samp} or (if 
#' "minnaert_full" is the chosen method) to \code{minnaert_samp}
#' @return The topographically corrected image.
#' @examples
#' # Mosaic the four ASTER DEM tiles needed to cover the Landsat image
#' ASTER_V002_LL <- raster(system.file('extdata/ASTER_V002_LL.dat', 
#' package='teamr'))
#' ASTER_V002_LR <- raster(system.file('extdata/ASTER_V002_LR.dat', 
#' package='teamr'))
#' ASTER_V002_UL <- raster(system.file('extdata/ASTER_V002_UL.dat', 
#' package='teamr'))
#' ASTER_V002_UR <- raster(system.file('extdata/ASTER_V002_UR.dat', 
#' package='teamr'))
#' DEM_mosaic <- mosaic(ASTER_V002_LL, ASTER_V002_LR, ASTER_V002_UR, 
#' ASTER_V002_UL, fun='mean')
#' 
#' # Crop and extend the DEM mosaic to match the Landsat image
#' L5TSR_1986_file <- system.file('extdata/L5TSR_1986.dat', package='teamr')
#' L5TSR_1986 <- stack(L5TSR_1986_file)
#' matched_DEM <- match_rasters(L5TSR_1986, DEM_mosaic)
#' 
#' # Read sun elevation and sun azimuth from .txt metadata file accompanying
#' # the Landsat file (as output from teampy Python package).
#' metadatafile <- extension(L5TSR_1986_file, 'txt')
#' sunelev <- 90 - as.numeric(get_metadata_item(metadatafile, 'SolarZenith'))
#' sunazimuth <- as.numeric(get_metadata_item(metadatafile, 'SolarAzimuth'))
#' 
#' slopeaspect <- slopeasp_par(matched_DEM)
#' # Apply the topographic correction
#' L5TSR_1986_topocorr <- topographic_corr(L5TSR_1986, sunelev, sunazimuth,
#'                                         slopeaspect, method='minslope')
#' 
#' plotRGB(L5TSR_1986, stretch='lin', r=3, g=2, b=1)
#' 
#' plotRGB(L5TSR_1986_topocorr, stretch='lin', r=3, g=2, b=1)
topographic_corr <- function(x, sunelev, sunazimuth, slopeaspect, method, 
                             filename=NULL, overwrite=FALSE, inparallel=FALSE, 
                             usesample=FALSE, ...) {
    if (!(class(x) %in% c('RasterLayer', 'RasterStack', 'RasterBrick'))) {
        stop('x must be a Raster* object')
    }
    if (!(class(slopeaspect) %in% c('RasterBrick', 'RasterStack'))) {
        stop('slopeaspect must be a RasterBrick or RasterStack object')
    }
    # Need to convert slope and aspect to SpatialGridDataFrame objects for 
    # topocorr. TODO: rewrite topocorr to handle RasterLayers
    slope <- as(raster(slopeaspect, layer=1), 'SpatialGridDataFrame')
    aspect <- as(raster(slopeaspect, layer=2), 'SpatialGridDataFrame')
    if (inparallel == TRUE) {
        # Set layer to NULL to pass R CMD CHECK without notes
        layer=NULL
        corr_img <- foreach(layer=unstack(x), .combine='addLayer', 
                            .multicombine=TRUE, .init=raster(), 
                            .packages=c('raster', 'rgdal', 'landsat')) %dopar% {
            img_df <- as(layer, 'SpatialGridDataFrame')
            if (method == 'minnaert_full') {
                minnaert_data <- minnaert_samp(img_df, slope, aspect, 
                                               sunelev=sunelev, 
                                               sunazimuth=sunazimuth, 
                                               usesample=usesample, ...)
                corr_df <- minnaert_data$minnaert
            } else {
                corr_df <- topocorr(img_df, slope, aspect, sunelev=sunelev, 
                                    sunazimuth=sunazimuth, method, usesample,
                                    ...)
            }
            raster(corr_df)
        }
    } else {
        corr_layers <- c()
        for (layer_num in 1:nlayers(x)) {
            message(paste0('Running topocorr on layer ', layer_num, ' of ', nlayers(x), '...'))
            img_df <- as(raster(x, layer=layer_num), 'SpatialGridDataFrame')
            if (method == 'minnaert_full') {
                minnaert_data <- minnaert_samp(img_df, slope, aspect, 
                                               sunelev=sunelev, 
                                               sunazimuth=sunazimuth, 
                                               usesample=usesample, ...)
                corr_df <- minnaert_data$minnaert
            } else {
                corr_df <- topocorr(img_df, slope, aspect, sunelev=sunelev, 
                                    sunazimuth=sunazimuth, method, ...)
            }
            corr_layers <- c(corr_layers, list(raster(corr_df)))
        }
        corr_img <- stack(corr_layers)
    }
    if (!is.null(filename)) {
        writeRaster(corr_img, filename, overwrite=overwrite)
    }
    return(corr_img)
}
