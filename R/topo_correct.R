
PLOT_WIDTH <- 3.25
PLOT_HEIGHT <- 3.25
DPI <- 300

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

mosaicDEMs <- function(DEMfilelist, ...) {
    require(raster)
    message(paste('Mosaicing DEMs with', DEMfilelist[1], 'as base file'))
    for (DEMfile in DEMfilelist[2:]) {
        # TODO: Need to ensure that projection systems match
        mosaic_img <- mosaic(mosaic_img, raster(DEMfile), fun=mean)
    }
    return(mosaic_img)
}

matchrasters <- function(base_image_file, match_image_file) {
    base_image <- raster(base_image_file)
    match_image <- raster(match_image_file)
}
base_image <- 'H:/Data/TEAM/VB/Rasters/Landsat/1986_037_LT5/proc/lndsr.LT50150531986037XXX18.bsq'
match_image <- 
matchrasters(base_image, match_image)

get_metadata_item <- function(metadatafile, item) {
    if (!file.exists(metadatafile) {
        stop(paste('Could not find metadata file', metadatafile))
    }
    metadata <- read.csv(metadatafile, stringsAsFactors=FALSE)
    rownum <- which(grepl(item, metadata$item))
    return(metadata$value[rownum])
}
#metadatafile <- 'H:/Data/TEAM/VB/Rasters/Landsat/1986_037_LT5/proc/lndsr.LT50150531986037XXX18.txt'
#get_metadata_item(metadatafile, 'SolarZenith')

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
