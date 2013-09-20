#' Makes a polygon shapefile of image extent
#'
#' @export
#' @param raster_file Path to the raster file.
#' @param out_file Filename for the output shapefile.
svm_classify <- function(raster_file, out_file) {
}

# shp <- 'H:/Data/TEAM/VB/Vectors/VB_training_1986_037_LT5.shp'
# img_file <- 'H:/Data/TEAM/VB/Rasters/Landsat/1986_037_LT5/proc/lndsr.LT50150531986037XXX18.bsq'

# require(teamr)
# shp <- 'H:/Data/TEAM/teamr_data/L5TSR_1986_training.shp'
# img_file <- 'H:/Data/TEAM/teamr_data/L5TSR_1986.dat'
# rast <- brick(img_file)
# rast <- stack(img_file, bands=as.integer(c(1, 2, 3, 4)))
# training_data <- extract_training_data(rast, shp)
# 
# testSVM <- svm(y ~ ., data=training_data)
# 
# tune.svm()
