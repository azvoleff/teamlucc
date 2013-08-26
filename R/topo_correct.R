PLOT_WIDTH <- 3.25
PLOT_HEIGHT <- 3.25
DPI <- 300

get_metadata_item <- function(metadatafile, item) {
    if (!file.exists(metadatafile)) {
        stop(paste('Could not find metadata file', metadatafile))
    }
    metadata <- read.csv(metadatafile, stringsAsFactors=FALSE)
    rownum <- which(grepl(item, metadata$item))
    return(metadata$value[rownum])
}
#metadatafile <- 'H:/Data/TEAM/VB/Rasters/Landsat/1986_037_LT5/proc/lndsr.LT50150531986037XXX18.txt'
#get_metadata_item(metadatafile, 'SolarZenith')

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

mosaicimgs <- function(baseimg, ...) {
    require(raster)
    require(rgdal)
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
library(raster)
DEM1 <- raster(system.file('extdata/ASTER_V002_LL.dat', package='LDPKR'))
DEM2 <- raster(system.file('extdata/ASTER_V002_LR.dat', package='LDPKR'))
DEM3 <- raster(system.file('extdata/ASTER_V002_UR.dat', package='LDPKR'))
DEM4 <- raster(system.file('extdata/ASTER_V002_UL.dat', package='LDPKR'))
DEM_mosaic_img <- mosaicimgs(DEM1, list(DEM2, DEM3, DEM4))

matchrasters <- function(baseimg, matchimg) {
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
    message('Resampling matchimg to base...')
    #resample(
    return(outimg)
}
library(raster)
DEM1 <- raster(system.file('extdata/ASTER_V002_LL.dat', package='LDPKR'))
DEM2 <- raster(system.file('extdata/ASTER_V002_LR.dat', package='LDPKR'))
DEM3 <- raster(system.file('extdata/ASTER_V002_UR.dat', package='LDPKR'))
DEM4 <- raster(system.file('extdata/ASTER_V002_UL.dat', package='LDPKR'))
DEM_mosaic_img <- mosaicimgs(DEM1, list(DEM2, DEM3, DEM4))

baseimg <- DEM_mosaic_img
matchimg <- raster(system.file('extdata/ASTER_V002_UL.dat', package='LDPKR'))
matched_image <- matchrasters(baseimg, matchimg)

topographic_corr <- function(imgfile, DEMfile, ...) {
    require(raster)
    require(landsat)
    img <- raster(imgfile)
    DEM <- raster(DEMfile)
    message('Calculating slope and aspect...')
    DEM.slopeasp <- slopeasp(DEM)
    sunelev <- as.numeric(get_metadata_item(extension(imgfile, 'txt'), 
                                            'SolarZenith'))
    sunazimuth <- as.numeric(get_metadata_item(extension(imgfile, 'txt'), 
                                               'SolarAzimuth'))
    message('Performing topographic correction...')
    corr_img <- topocorr(img, DEM.slopeasp$slope, DEM.slopeasp$aspect, 
                         sunelev=sunelev, sunazimuth=sunazimuth, ...)
    return(corr_img)
}
