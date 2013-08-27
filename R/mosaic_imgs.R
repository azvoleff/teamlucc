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
