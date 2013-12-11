#' Calculate change-trajectory image
#'
#' @export
#' @import spatial.tools
#' @param initial initial cover class as \code{RasterLayer}
#' @param chg_mag change magnitude \code{RasterLayer} from \code{CVAPS}
#' @param chg_dir change direction \code{RasterLayer} from \code{CVAPS}
#' @param threshold the threshold to use as a minimum when determining change 
#' areas (can use \code{DFPS} to determine this value).
#' @param filename filename to save the output \code{RasterLayer} to disk 
#' (optional)
#' @param classnames an optional vector of classnames to output with the 
#' returned trajectory lookup table
#' @param ignorepersistence whether to ignore persistence of a class (if 
#' ignored all pixels where a class persists will be set to NA)
#' @return a table of all possible trajectories, with their \code{classnames} 
#' (if specified) and the integer codes used to indicate specific trajectories 
#' in the output image.
#' @references Chen, J., P. Gong, C.  He, R.  Pu, and P.  Shi.  2003.
#' Land-use/land-cover change detection using improved change-vector analysis.
#' Photogrammetric Engineering and Remote Sensing 69:369-380.
#' 
#' Chen, J., X. Chen, X. Cui, and J. Chen. 2011. Change vector analysis in
#' posterior probability space: a new method for land cover change detection.
#' IEEE Geoscience and Remote Sensing Letters 8:317-321.
chg_traj <- function(initial, chg_mag, chg_dir, threshold, filename=NULL, 
                     classnames=NULL, ignorepersistence=TRUE) {
    if (proj4string(initial) != proj4string(chg_mag) ) {
        stop('initial and chg_mag coordinate systems do not match')
    } else if (proj4string(initial) != proj4string(chg_dir) ) {
        stop('initial and chg_dir coordinate systems do not match')
    }
    if (extent(initial) != extent(chg_mag)) {
        stop('extent of initial does not match extent of chg_mag')
    } else if (extent(initial) != extent(chg_dir)) {
        stop('extent of initial does not match extent of chg_dir')
    }
    if (nlayers(initial) > 1) stop('initial has more than 1 layer')
    if (nlayers(chg_mag) > 1) stop('chg_mag has more than 1 layer')
    if (nlayers(chg_dir) > 1) stop('chg_dir has more than 1 layer')

    # Make a lookup table of codes for each type of transition
    classcodes <- sort(unique(getValues(initial)))
    traj_lut <- expand.grid(t0_code=classcodes, t1_code=classcodes)
    if (!is.null(classnames)) {
        if (length(classnames) != length(classcodes)) {
            stop('classnames must be NULL or a vector of length equal to number of classes in initial image')
        }
        traj_lut$t0_name <- classnames[match(traj_lut$t0_code, classcodes)]
        traj_lut$t1_name <- classnames[match(traj_lut$t1_code, classcodes)]
    }
    if (ignorepersistence) traj_lut <- traj_lut[!(traj_lut$t0_code == traj_lut$t1_code),]
    # Code trajectories by summing t0 and t1 after multiplying t1 by the number 
    # of classes.
    traj_lut$Code <- traj_lut$t0_code + traj_lut$t1_code * length(classcodes)

    calc_chg_traj <- function(initial, chg_mag, chg_dir, classcodes, threshold, 
                              ignorepersistence, ...) {
        # Code trajectories by summing t0 and chg_dir after multiplying chg_dir 
        # by the number of classes.
        traj <- initial + chg_dir * length(classcodes)
        traj[chg_mag < threshold] <- NA
        # Ignore persistence of a class if desired
        if (ignorepersistence) traj[initial == chg_dir] <- NA
        traj <- array(traj, dim=c(dim(initial)[1], dim(initial)[2], 1))
        return(traj)
    }
    out <- rasterEngine(initial=initial, chg_mag=chg_mag, chg_dir=chg_dir, 
                        fun=calc_chg_traj, args=list(classcodes=classcodes, 
                                                     threshold=threshold,
                                                     ignorepersistence=ignorepersistence), 
                        filename=filename)

    return(list(traj_lut=traj_lut, chg_traj=out))
}
