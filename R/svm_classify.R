#' Makes a polygon shapefile of image extent
#'
#' @export
#' @param raster_file Path to the raster file.
#' @param out_file Filename for the output shapefile.
svm_classify <- function(raster_file, out_file) {
}

# shp <- 'H:/Data/TEAM/VB/Vectors/VB_training_1986_037_LT5.shp'
# img_file <- 'H:/Data/TEAM/VB/Rasters/Landsat/1986_037_LT5/proc/lndsr.LT50150531986037XXX18.bsq'

require(teamr)
t0_shp <- 'H:/Data/TEAM/teamr_data/L5TSR_1986_training.shp'
t0_img <- stack('H:/Data/TEAM/teamr_data/L5TSR_1986.dat',
                bands=as.integer(c(1, 2, 3, 4)))
t1_shp <- 'H:/Data/TEAM/teamr_data/L5TSR_2001_training.shp'
t1_img <- stack('H:/Data/TEAM/teamr_data/L5TSR_2001.dat',
                bands=as.integer(c(1, 2, 3, 4)))

svm_t0 <- tune.svm(y ~ ., data=extract_training_data(t0_img, t0_shp))
svm_t1 <- tune.svm(y ~ ., data=extract_training_data(t1_img, t1_shp))

svm_t0
svm_t1

tune.svm()
