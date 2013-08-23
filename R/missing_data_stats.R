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

image_list$pct_cloud <- NA
image_list$pct_missing <- NA

for (n in 1:nrow(image_list)) {
    in_prefix <- sub('.hdf$', '', image_list[n, ]$file_path)
    out_prefix <- sub('orig', 'proc', in_prefix)
    print(paste('Processing', in_prefix))

    fill_QA_out_prefix <- paste(out_prefix, '_fill_QA_crop', sep='')
    fill_QA <- raster(paste(in_prefix, '_fill_QA.bsq', sep=''))
    fill_QA_masked <- mask_ZOI(fill_QA, zoi)
    fill_QA_masked_stats <- table(getValues(fill_QA_masked))
    image_list$pct_missing[n] <- 100 - as.numeric(fill_QA_masked_stats['0'] / 
                                                  sum(fill_QA_masked_stats) * 
                                                  100)
     
    cloud_out_prefix <- paste(out_prefix, '_cloud_crop', sep='')
    cloud <- raster(paste(out_prefix, '_comb_cloud_mask.bsq', sep=''))
    cloud_crop_masked <- mask_ZOI(cloud, zoi)
    cloud_crop_masked_stats <- table(getValues(cloud_crop_masked))
    image_list$pct_cloud[n] <- as.numeric(cloud_crop_masked_stats['255'] / 
                                          sum(cloud_crop_masked_stats) * 100)
}

write.csv(image_list, file='H:/Data/TEAM/VB/Rasters/Landsat/image_list.csv', 
          row.names=FALSE)
