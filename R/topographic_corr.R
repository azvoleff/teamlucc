#' Topographically correct a raster
#'
#' Performs topographic correction using code based on \code{\link{topocorr}} 
#' from the \code{Landsat} package. The code in this package has been modifed 
#' from \code{\link{topocorr}} to allow using a subsample of the image for 
#' Minnaert k calculations, and to provide the option of running the 
#' topographic correction in parallel using \code{\link{foreach}}.
#'
#' @export
#' @import foreach
#' @param x an image to correct
#' @param slopeaspect a \code{RasterBrick} or \code{RasterStack} with two 
#' layers.  First layer should be the slope, second layer should be aspect. The 
#' slope and aspect are defined as in \code{slopeasp} in the \code{landsat} 
#' package.  \code{\link{slopeasp_par}} will output the slope and aspect using 
#' the proper definition and as a \code{RasterBrick}.
#' @param sunelev sun elevation in degrees
#' @param sunazimuth sun azimuth in degrees
#' @param method the topographic correction method to use. See the help for 
#' \code{\link{topocorr}}.
#' @param filename file on disk to save \code{Raster*} to (optional)
#' @param overwrite whether to overwrite \code{filename} if it already exists
#' @param inparallel whether to run correction in parallel using \code{foreach}
#' @param sampleindices (optional) row-major indices of sample pixels to use in 
#' regression models used for some topographic correction methods (like 
#' Minnaert). Useful when handling very large images. See
#' \code{\link{gridsample}} for one method of calculating these indices.
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
topographic_corr <- function(x, slopeaspect, sunelev, sunazimuth, method, 
                             filename=NULL, overwrite=FALSE, inparallel=FALSE, 
                             sampleindices=NULL) {
    if (!(class(x) %in% c('RasterLayer', 'RasterStack', 'RasterBrick'))) {
        stop('x must be a Raster* object')
    }
    if (!(class(slopeaspect) %in% c('RasterBrick', 'RasterStack'))) {
        stop('slopeaspect must be a RasterBrick or RasterStack object')
    }
    # Need to convert slope and aspect to SpatialGridDataFrame objects for 
    # topocorr. TODO: rewrite topocorr to handle RasterLayers
    slope <- raster(slopeaspect, layer=1)
    aspect <- raster(slopeaspect, layer=2)
    if (inparallel == TRUE && (nlayers(x) > 1)) {
        # Set layer to NULL to pass R CMD CHECK without notes
        uncorr_layer=NULL
        corr_img <- foreach(uncorr_layer=unstack(x), .combine='addLayer', 
                            .multicombine=TRUE, .init=raster(), 
                            .packages=c('raster', 'rgdal', 'teamr')) %dopar% {
            if (method == 'minnaert_full') {
                minnaert_data <- minnaert_samp(uncorr_layer, slope, aspect, 
                                               sunelev=sunelev, 
                                               sunazimuth=sunazimuth, 
                                               sampleindices=sampleindices)
                corr_layer <- minnaert_data$minnaert
            } else {
                corr_layer <- topocorr_samp(uncorr_layer, slope, aspect, 
                                            sunelev=sunelev, sunazimuth=sunazimuth, 
                                            method=method, sampleindices=sampleindices)
            }
        }
    } else {
        corr_layers <- c()
        for (layer_num in 1:nlayers(x)) {
            message(paste0('Running topocorr on layer ', layer_num, ' of ', nlayers(x), '...'))
            if (nlayers(x) > 1) {
                uncorr_layer <- raster(x, layer=layer_num)
            } else {
                uncorr_layer <- x
            }
            if (method == 'minnaert_full') {
                minnaert_data <- minnaert_samp(uncorr_layer, slope, aspect, 
                                               sunelev=sunelev, 
                                               sunazimuth=sunazimuth, 
                                               sampleindices=sampleindices)
                corr_layer <- minnaert_data$minnaert
            } else {
                corr_layer <- topocorr_samp(uncorr_layer, slope, aspect, 
                                            sunelev=sunelev, 
                                            sunazimuth=sunazimuth, 
                                            method=method, 
                                            sampleindices=sampleindices)
            }
            corr_layers <- c(corr_layers, list(corr_layer))
        }
        if (nlayers(x) > 1) {
            corr_img <- stack(corr_layers)
        } else {
            corr_img <- corr_layers[[1]]
        }
    }
    if (!is.null(filename)) {
        writeRaster(corr_img, filename, overwrite=overwrite)
    }
    return(corr_img)
}
