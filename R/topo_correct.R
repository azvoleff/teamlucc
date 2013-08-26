get_metadata_item <- function(metadatafile, item) {
    if (!file.exists(metadatafile)) {
        stop(paste('Could not find metadata file', metadatafile))
    }
    metadata <- read.csv(metadatafile, stringsAsFactors=FALSE)
    rownum <- which(grepl(item, metadata$item))
    return(metadata$value[rownum])
}

get_band_names_from_hdr <- function(hdr_file) {
    require(raster)
    txt <- readLines(hdr_file)
    line_num <- which(grepl('^band names', a)) + 1
    band_names <- c()
    for (n in line_num:length(txt)) {
        band_names <- c(band_names, gsub('[,}][[:space:]]*$', '', txt[n]))
        if (grepl('}', txt[n])) {
            break
        }
    }
    return(band_names)
}

mosaic_imgs <- function(baseimg, ...) {
    require(raster)
    dots <- list(...)
    if (length(dots) < 1) {
        stop('no mosaic images provided')
    } else if (is.list(dots[[1]])) {
        img_list <- dots[[1]]
        if (length(dots) > 1) {
            warning('second argument is a list - but only mosaicing two images')
        }
    } else {
        img_list <- dots
    }
    mosaic_img <- baseimg
    for (img in img_list) {
        message(paste('Adding image to mosaic...'))
        # TODO: Need to ensure that projection systems match
        mosaic_img <- mosaic(mosaic_img, img, fun=mean)
    }
    return(mosaic_img)
}
# library(raster)
# DEM1 <- raster(system.file('extdata/ASTER_V002_LL.dat', package='LDPKR'))
# DEM2 <- raster(system.file('extdata/ASTER_V002_LR.dat', package='LDPKR'))
# DEM3 <- raster(system.file('extdata/ASTER_V002_UR.dat', package='LDPKR'))
# DEM4 <- raster(system.file('extdata/ASTER_V002_UL.dat', package='LDPKR'))
# DEM_mosaic_img <- mosaic_imgs(DEM1, list(DEM2, DEM3, DEM4))

match_rasters <- function(baseimg, matchimg) {
    require(raster)
    if (projection(baseimg) != projection(matchimg)) {
        message('Coordinate systems do not match - reprojecting matchimg...')
        matchimg <- projectRaster(matchimg, baseimg)
    }
    # First crop out any overlapping area
    message('Cropping matchimg to base...')
    outimg <- crop(matchimg, baseimg)
    # Now extend borders of cropped raster to match base raster
    message('Extending matchimg to base...')
    outimg <- extend(matchimg, outimg)
    #message('Resampling matchimg to base...')
    #resample(outimg, baseimg)
    return(outimg)
}
# library(raster)
# # Mosaic the four ASTER DEM tiles needed to cover the Landsat image
# DEM1 <- raster(system.file('extdata/ASTER_V002_LL.dat', package='LDPKR'))
# DEM2 <- raster(system.file('extdata/ASTER_V002_LR.dat', package='LDPKR'))
# DEM3 <- raster(system.file('extdata/ASTER_V002_UR.dat', package='LDPKR'))
# DEM4 <- raster(system.file('extdata/ASTER_V002_UL.dat', package='LDPKR'))
# DEM_mosaic <- mosaic_imgs(DEM1, list(DEM2, DEM3, DEM4))
# 
# # Crop and extend the DEM mosaic to match the Landsat image
# L5TSR_1986 <- raster(system.file('extdata/L5TSR_1986.dat', package='LDPKR'))
# matched_DEM <- match_rasters(L5TSR_1986, DEM_mosaic)

topographic_corr <- function(img, DEM, sunelev, sunazimuth, ...) {
    require(raster)
    require(landsat)
    if (class(img) == 'SpatialGridDataFrame') {
        stop('img must be a Raster* object')
    }
    if (class(DEM) == 'SpatialGridDataFrame') {
        stop('DEM must be a Raster* object')
    } else {
        DEM_df <- as(DEM, "SpatialGridDataFrame")
    }
    message('Calculating slope and aspect...')
    DEM_slopeasp <- slopeasp(DEM_df)
    corr_img <- raster()
    for (bandnum in 1:nlayers(img)) {
        message(paste('Performing topographic correction on band ', bandnum, 
                      '...', sep=''))
        img_df <- as(raster(img, layer=bandnum), 'SpatialGridDataFrame')
        corr_df <- topocorr(img_df, DEM_slopeasp$slope, DEM_slopeasp$aspect, 
                             sunelev=sunelev, sunazimuth=sunazimuth, ...)
        corr_img <- addLayer(corr_img, raster(corr_df))
    }
    return(corr_img)
}
# library(raster)
# # Mosaic the four ASTER DEM tiles needed to cover the Landsat image
# DEM1 <- raster(system.file('extdata/ASTER_V002_LL.dat', package='LDPKR'))
# DEM2 <- raster(system.file('extdata/ASTER_V002_LR.dat', package='LDPKR'))
# DEM3 <- raster(system.file('extdata/ASTER_V002_UR.dat', package='LDPKR'))
# DEM4 <- raster(system.file('extdata/ASTER_V002_UL.dat', package='LDPKR'))
# DEM_mosaic <- mosaic_imgs(DEM1, list(DEM2, DEM3, DEM4))
# 
# # Crop and extend the DEM mosaic to match the Landsat image
# L5TSR_1986_file <- system.file('extdata/L5TSR_1986.dat', package='LDPKR')
# L5TSR_1986 <- stack(L5TSR_1986_file)
# matched_DEM <- match_rasters(L5TSR_1986, DEM_mosaic)
# 
# # Read sun elevation and sun azimuth from .txt metadata file accompanying the 
# # Landsat file (as output from LDPK Python tools)
# metadatafile <- extension(L5TSR_1986_file, 'txt')
# sunelev <- as.numeric(get_metadata_item(metadatafile, 'SolarZenith'))
# sunazimuth <- as.numeric(get_metadata_item(metadatafile, 'SolarAzimuth'))
# 
# # Apply the topographic correction
# L5TSR_1986_topocorr<- topographic_corr(L5TSR_1986, matched_DEM, sunelev,
#                                         sunazimuth, method='minslope')
# 
# plotRGB(L5TSR_1986, stretch='lin', r=3, g=2, b=1)
# 
# plotRGB(L5TSR_1986_topocorr, stretch='lin', r=3, g=2, b=1)
