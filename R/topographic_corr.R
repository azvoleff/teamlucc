#' Topographically correct a raster
#'
#' @export
#' @import landsat foreach
#' @param x An image to correct.
#' @param DEM A digital elevation model with the same coordinate system, 
#' extent, and resolution as /code{x}
#' @param sunazimuth Sun azimuth in degrees.
#' @param sunelev Sun elevation in degrees.
#' @param ... Additional arguments to be passed to \code{topocorr} in the 
#' \code{landsat} package.
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
#' # Apply the topographic correction
#' L5TSR_1986_topocorr <- topographic_corr(L5TSR_1986, matched_DEM, sunelev,
#'                                         sunazimuth, method='minslope')
#' 
#' plotRGB(L5TSR_1986, stretch='lin', r=3, g=2, b=1)
#' 
#' plotRGB(L5TSR_1986_topocorr, stretch='lin', r=3, g=2, b=1)
topographic_corr <- function(x, DEM, sunelev, sunazimuth, ...) {
    if (class(x) == 'SpatialGridDataFrame') {
        stop('x must be a Raster* object')
    }
    if (class(DEM) == 'SpatialGridDataFrame') {
        stop('DEM must be a Raster* object')
    } else {
        DEM_df <- as(DEM, "SpatialGridDataFrame")
    }
    message('Calculating slope and aspect...')
    DEM_slopeasp <- slopeasp_par(DEM)
    # Need to convert slope and aspect to SpatialGridDataFrame objects for 
    # topocorr. TODO: rewrite topocorr to handle RasterLayers
    slope <- as(DEM_slopeasp$slope, 'SpatialGridDataFrame')
    aspect <- as(DEM_slopeasp$aspect, 'SpatialGridDataFrame')
    message('Performing topographic correction...')
    corr_img <- foreach(layer=unstack(x), .combine='addLayer', .multicombine=TRUE, 
             .init=raster(), .packages=c('raster', 'rgdal', 'landsat')) %dopar% {
        img_df <- as(layer, 'SpatialGridDataFrame')
        corr_df <- topocorr(img_df, slope, aspect, 
                             sunelev=sunelev, sunazimuth=sunazimuth, ...)
        raster(corr_df)
    }
    return(corr_img)
}
