#' Topographically correct a raster
#'
#' @export
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
#' data(ASTER_V002_LL)
#' data(ASTER_V002_LR)
#' data(ASTER_V002_UR)
#' data(ASTER_V002_UL)
#' DEM_mosaic <- mosaic(ASTER_V002_LL, ASTER_V002_LR, ASTER_V002_UR, ASTER_V002_UL, fun='mean')
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
#' L5TSR_1986_topocorr<- topographic_corr(L5TSR_1986, matched_DEM, sunelev,
#'                                         sunazimuth, method='minslope')
#' 
#' plotRGB(L5TSR_1986, stretch='lin', r=3, g=2, b=1)
#' 
#' plotRGB(L5TSR_1986_topocorr, stretch='lin', r=3, g=2, b=1)
topographic_corr <- function(x, DEM, sunelev, sunazimuth, ...) {
    require(raster)
    require(landsat)
    if (class(x) == 'SpatialGridDataFrame') {
        stop('x must be a Raster* object')
    }
    if (class(DEM) == 'SpatialGridDataFrame') {
        stop('DEM must be a Raster* object')
    } else {
        DEM_df <- as(DEM, "SpatialGridDataFrame")
    }
    message('Calculating slope and aspect...')
    DEM_slopeasp <- slopeasp(DEM_df)
    corr_img <- raster()
    for (bandnum in 1:nlayers(x)) {
        message(paste('Performing topographic correction on band ', bandnum, 
                      '...', sep=''))
        img_df <- as(raster(x, layer=bandnum), 'SpatialGridDataFrame')
        corr_df <- topocorr(img_df, DEM_slopeasp$slope, DEM_slopeasp$aspect, 
                             sunelev=sunelev, sunazimuth=sunazimuth, ...)
        corr_img <- addLayer(corr_img, raster(corr_df))
    }
    return(corr_img)
}
