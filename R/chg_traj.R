#' Calculate change-trajectory image
#'
#' @export
#' @importFrom spatial.tools rasterEngine
#' @param initial initial cover class as \code{RasterLayer}
#' @param chg_mag change magnitude \code{RasterLayer} from \code{CVAPS}
#' @param chg_dir change direction \code{RasterLayer} from \code{CVAPS}
#' @param chg_threshold the threshold to use as a minimum when determining change 
#' areas (can use \code{DFPS} to determine this value).
#' @param filename filename to save the output \code{RasterLayer} to disk 
#' (optional)
#' @param classnames an optional vector of classnames to output with the 
#' returned trajectory lookup table
#' @param overwrite whether to overwrite existing files (otherwise an error 
#' will be raised)
#' @param overwrite whether to overwrite existing files (otherwise an error 
#' will be raised)
#' @param ... additional parameters to pass to rasterEngine
#' @return a table of all possible trajectories, with their \code{classnames} 
#' (if specified) and the integer codes used to indicate specific trajectories 
#' in the output image.
#' @details Processing can be done in parallel using all using the cluster 
#' facilities in the \code{spatial.tools} package. To enable clustering, call 
#' \code{beginCluster} before running \code{classify}.  To stop the 
#' cluster when finished, call \code{endCluster}.
#' @references Chen, J., P. Gong, C.  He, R.  Pu, and P.  Shi.  2003.
#' Land-use/land-cover change detection using improved change-vector analysis.
#' Photogrammetric Engineering and Remote Sensing 69:369-380.
#' 
#' Chen, J., X. Chen, X. Cui, and J. Chen. 2011. Change vector analysis in
#' posterior probability space: a new method for land cover change detection.
#' IEEE Geoscience and Remote Sensing Letters 8:317-321.
chg_traj <- function(initial, chg_mag, chg_dir, chg_threshold, filename,
                     classnames=NULL, overwrite=FALSE, ...) {
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
    if (!missing(filename) && file_test('-f', filename) && !overwrite) {
        stop('output file already exists and overwrite=FALSE')
    }

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
    # Code trajectories by summing t0 and t1 after multiplying t1 by the number 
    # of classes.
    traj_lut$Code <- traj_lut$t0_code + traj_lut$t1_code * length(classcodes)

    calc_chg_traj <- function(initial, chg_mag, chg_dir, classcodes, 
                              chg_threshold, ...) {
        # Code change trajectories by summing t0 and chg_dir after multiplying 
        # chg_dir by the number of classes
        traj <- initial + chg_dir * length(classcodes)
        # Code classes that persist, being sure to ignore NAs
        persist_pixels <- which(chg_mag < chg_threshold)
        persist_pixels <- persist_pixels[!is.na(persist_pixels)]
        traj[persist_pixels] <- initial[persist_pixels] + initial[persist_pixels] * length(classcodes)
        traj <- array(traj, dim=c(dim(initial)[1], dim(initial)[2], 1))
        return(traj)
    }
    out <- rasterEngine(initial=initial, chg_mag=chg_mag, chg_dir=chg_dir, 
                        fun=calc_chg_traj,
                        args=list(classcodes=classcodes, 
                                  chg_threshold=chg_threshold),
                        datatype='INT2S', ...)

    # spatial.tools can only output the raster package grid format - so output 
    # to a tempfile in that format then copy over to the requested final output 
    # format if a filename was supplied
    if (!missing(filename)) {
        out <- writeRaster(out, filename=filename, overwrite=overwrite, 
                           datatype='INT2S')
    }

    return(list(lut=traj_lut, traj=out))
}
