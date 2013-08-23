library(rgdal)
library(raster)
library(rgeos) # needed for gBuffer function

image_list <- read.csv('H:/Data/TEAM/VB/Rasters/Landsat/image_list.csv')
zoi <- readOGR('H:/Data/TEAM/VB/Shapefiles', 'VB_ZOI_GEO')

mask_ZOI <- function(img, zoi, bufferwidth=1000) {
    zoi <- spTransform(zoi, CRS(proj4string(img)))
    # Buffer the ZOI so that cloud fill algorithms can have room to work when 
    # interpolating missing data
    zoi_buffer <- gBuffer(zoi, width=bufferwidth)
    img_crop <- crop(img, zoi_buffer)
    img_crop_mask <- mask(img_crop, zoi_buffer)
}

band_names <- list('band1', 'band2', 'band3', 'band4', 'band5', 'band6', 
                   'band7', 'comb_cloud_mask', 'fill_QA')

for (n in 1:nrow(image_list)) {
    in_prefix <- sub('.hdf$', '', image_list[n, ]$file_path)
    out_prefix <- sub('orig', 'proc', in_prefix)
    print(paste('Processing', in_prefix))

    for (band in band_names) {
        b_out_prefix <- paste(out_prefix, band, 'zoi_crop', sep='_')
        if (file.exists(paste(b_out_prefix, '.envi', sep=''))) {
            print(paste('skipping', band, '- file already exists'))
            next
        }
        print(paste('writing', band))
        b <- raster(paste(in_prefix, '_', band, '.bsq', sep=''))
        b_masked <- mask_ZOI(b, zoi)
        writeRaster(b_masked, b_out_prefix, format='ENVI')
    }

}
