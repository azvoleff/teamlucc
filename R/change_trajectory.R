#' Calculate change-trajectory image
#'
#' @export
#' @param initial_img initial cover class as \code{RasterLayer}
#' @param chg_img change direction \code{RasterLayer} from \code{CVAPS}
#' @param out_file_base the base name to use when naming the output files
#' @param threshold the threshold to use as a minimum when determining change 
#' areas (can use \code{DFPS} to determine this value).
#' @return a \code{RasterLayer} where no change areas are NA, and change areas 
#' are coded as
#' @references Chen, J., P. Gong, C. He, R. Pu, and P. Shi. 2003.
#' Land-use/land-cover change detection using improved change-vector analysis.
#' Photogrammetric Engineering and Remote Sensing 69:369-380.
#' 
#' Chen, J., X. Chen, X. Cui, and J. Chen. 2011. Change vector analysis in
#' posterior probability space: a new method for land cover change detection.
#' IEEE Geoscience and Remote Sensing Letters 8:317-321.
change_trajectory <- function(initial_img, chg_img, out_file_base, threshold) {
    if (proj4string(initial_img) != proj4string(chg_img)) {
        stop('Error: t0 and t1 coordinate systems do not match')
    }
    if (extent(initial_img) != extent(chg_img)) {
        stop('Error: t0 and t1 extents do not match')
    }
    if ((nlayers(initial_img) > 1) | (nlayers(initial_img)) > 1) {
        stop('Error: neither initial_img nor chg_img should have more than 1 layer')
    }
    # Make a lookup table of codes for each type of transition
    classes <- unique(getValues(initial_img))
    traj_lut <- expand.grid(t0=classes, t1=classes)
    traj_lut$Code <- seq(1, nrow(traj_lut))
    out_traj <- raster(initial_img)
    out_traj <- writeStart(out_traj, paste(out_file_base, '_chgtraj.envi', sep=''))
    bs <- blockSize(initial_img)
    for (block_num in 1:bs$n) {
        initial_img_block <- getValues(initial_img, row=bs$row[block_num], 
                                       nrows=bs$nrows[block_num])
        chg_img_block <- getValues(chg_img, row=bs$row[block_num], 
                                   nrows=bs$nrows[block_num])
        traj <- initial_img_block
        traj[traj < threshold] <- NA
        #TODO: Finish this code
        #chg_pixels < [traj >= threshold]
        for (t0_n in 1:nrow(traj_lut)) {
            for (t1_n in 1:nrow(traj_lut)) {
                #TODO: Finish this code
                traj[initial_img_block > threshold]
            }
        }
                                   
        traj[chg_img > threshold]
    }
    out_traj <- writeStop(out_traj)
}
