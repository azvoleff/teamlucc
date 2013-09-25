#' Calculate change-trajectory image
#'
#' @export
#' @param initial initial cover class as \code{RasterLayer}
#' @param chg_dir change direction \code{RasterLayer} from \code{CVAPS}
#' @param chg_mag change magnitude \code{RasterLayer} from \code{CVAPS}
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
change_trajectory <- function(initial, chg_mag, chg_dir, out_file_base, 
                              threshold, classnames=NULL, 
                              ignorepersistence=FALSE) {
    if (proj4string(initial) != proj4string(chg_mag) ) {
        stop('Error: initial and chg_mag coordinate systems do not match')
    } else if (proj4string(initial) != proj4string(chg_dir) ) {
        stop('Error: initial and chg_dir coordinate systems do not match')
    }
    if (extent(initial) != extent(chg_mag)) {
        stop('Error: extent of initial does not match extent of chg_mag')
    } else if (extent(initial) != extent(chg_dir)) {
        stop('Error: extent of initial does not match extent of chg_dir')
    }
    if (nlayers(initial) > 1) stop('Error: initial has more than 1 layer')
    if (nlayers(chg_mag) > 1) stop('Error: chg_mag has more than 1 layer')
    if (nlayers(chg_dir) > 1) stop('Error: chg_dir has more than 1 layer')

    # Make a lookup table of codes for each type of transition
    classcodes <- sort(unique(getValues(initial)))
    traj_lut <- expand.grid(t0_code=classcodes, t1_code=classcodes)
    if (!is.null(classnames)) {
        if (length(classnames) != length(classcodes)) {
            stop('Error: classnames must be NULL or a vector of length equal to number of classes in initial image')
        }
        traj_lut$t0_name <- classnames[match(traj_lut$t0_code, classcodes)]
        traj_lut$t1_name <- classnames[match(traj_lut$t1_code, classcodes)]
    }
    if (ignorepersistence) traj_lut <- traj_lut[!(traj_lut$t0_code == traj_lut$t1_code),]
    # Code trajectories by summing t0 and t1 after multiplying t1 by the number 
    # of classes.
    traj_lut$Full_Traj_Code <- traj_lut$t0_code + traj_lut$t1_code * length(classcodes)

    out_traj <- raster(initial)
    out_traj <- writeStart(out_traj, paste(out_file_base, '_chgtraj.envi', sep=''))
    bs <- blockSize(initial)
    for (block_num in 1:bs$n) {
        initial_blk <- getValues(initial, row=bs$row[block_num], 
                             nrows=bs$nrows[block_num])
        chg_mag_blk <- getValues(chg_mag, row=bs$row[block_num], 
                                 nrows=bs$nrows[block_num])
        chg_dir_blk <- getValues(chg_dir, row=bs$row[block_num], 
                                 nrows=bs$nrows[block_num])
        # Code trajectories by summing t0 and t1 after multiplying t1 by the 
        # number of classes.
        traj <- initial_blk + chg_dir_blk * length(classcodes)
        traj[chg_mag_blk < threshold] <- NA
        # Ignore persistence of a class if desired
        if (ignorepersistence) traj[initial_blk == chg_dir_blk] <- NA
        writeValues(out_traj, traj, bs$row[block_num])
    }
    out_traj <- writeStop(out_traj)

    return(traj_lut)
}
