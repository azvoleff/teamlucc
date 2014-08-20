#' Calculate change-trajectory lookup table
#'
#' This function will format a lookup table (lut) to allow coding change 
#' trajectories. Useful for use in conjunction with \code{\link{chg_traj}}.
#'
#' @export
#' @param class_codes a list of integer codes used to code land use/cover 
#' classes
#' @param class_names an (optional) list of class names as character vectors
#' @examples
#' lut <- traj_lut(c(1, 2), c("NonForest", "Forest"))
traj_lut <- function(class_codes, class_names=NULL) {
    lut <- expand.grid(t0_code=class_codes, t1_code=class_codes)
    if (!is.null(class_names)) {
        if (length(class_names) != length(class_codes)) {
            stop('class_names must be NULL or a vector of length equal to number of classes in initial image')
        }
        lut$t0_name <- class_names[match(lut$t0_code, class_codes)]
        lut$t1_name <- class_names[match(lut$t1_code, class_codes)]
    }
    # Code trajectories by summing t0 and t1 after multiplying t1 by the number 
    # of classes.
    lut$Code <- lut$t0_code + lut$t1_code * length(class_codes)
    # Exclude classes that are persistence - CVAPS doesn't directly code the 
    # class for classes that persist
    lut <- lut[!(lut$t0_code == lut$t1_code), ]
    return(lut)
}

#' Calculate change-trajectory image
#'
#' This function will calculate trajectories of land cover change using the 
#' Change Vector Analysis in Posterior Probability Space (CVAPS) approach of 
#' comparing posterior probabilities of class membership with an automatically 
#' determined threshold. Areas of no change are coded as -1. A lookup table for 
#' the codes output by \code{chg_traj} can be calculated with \code{traj_lut}.
#'
#' This function will run in parallel if a parallel backend is registered with 
#' \code{\link{foreach}}.
#'
#' @export
#' @importFrom spatial.tools rasterEngine
#' @param chg_mag change magnitude \code{RasterLayer} from \code{CVAPS}
#' @param chg_dir change direction \code{RasterLayer} from \code{CVAPS}
#' @param chg_threshold the threshold to use as a minimum when determining change 
#' areas (can use \code{DFPS} to determine this value).
#' @param filename filename to save the output \code{RasterLayer} to disk 
#' (optional)
#' @param overwrite whether to overwrite existing files (otherwise an error 
#' will be raised)
#' @param ... additional parameters to pass to rasterEngine
#' @return a {RasterLayer} of change trajectories, with change trajectories 
#' coded as in the \code{lut} output by \code{traj_lut}
#' @references Chen, J., P. Gong, C.  He, R.  Pu, and P.  Shi.  2003.
#' Land-use/land-cover change detection using improved change-vector analysis.
#' Photogrammetric Engineering and Remote Sensing 69:369-380.
#' 
#' Chen, J., X. Chen, X. Cui, and J. Chen. 2011. Change vector analysis in
#' posterior probability space: a new method for land cover change detection.
#' IEEE Geoscience and Remote Sensing Letters 8:317-321.
#' @examples
#' \dontrun{
#' t0_train_data <- get_pixels(L5TSR_1986, L5TSR_1986_2001_training, "class_1986",training=.6)
#' t0_model <- train_classifier(t0_train_data)
#' t0_preds <- classify(L5TSR_1986, t0_model)
#' t1_train_data <- get_pixels(L5TSR_2001, L5TSR_1986_2001_training, "class_2001", training=.6)
#' t1_model <- train_classifier(t1_train_data)
#' t1_preds <- classify(L5TSR_2001, t1_model)
#' t0_t1_chgmag <- chg_mag(t0_preds$probs, t1_preds$probs)
#' t0_t1_chgdir <- chg_dir(t0_preds$probs, t1_preds$probs)
#' 
#' lut <- traj_lut(t0_preds$codes$code, t0_preds$codes$class)
#' t0_t1_chgtraj <- chg_traj(lut, t0_t1_chgmag, t0_t1_chgdir, .5)
#' 
#' # Change areas are coded following the above lookup-table (lut):
#' plot(t0_t1_chgtraj)
#' 
#' # No change areas are -1:
#' plot(t0_t1_chgtraj == -1)
#' }
chg_traj <- function(chg_mag, chg_dir, chg_threshold, filename, 
                     overwrite=FALSE, ...) {
    if (nlayers(chg_mag) > 1) stop('chg_mag has more than 1 layer')
    if (nlayers(chg_dir) > 1) stop('chg_dir has more than 1 layer')
    compareRaster(chg_mag, chg_dir)
    if (!missing(filename) && file_test('-f', filename) && !overwrite) {
        stop('output file already exists and overwrite=FALSE')
    }

    calc_chg_traj <- function(chg_mag, chg_dir, chg_threshold, ...) {
        # Trajectories in chg_dir were coded by summing t0 and t1 classes after 
        # multiplying t1 class by the number of classes
        chg_dir[chg_mag < chg_threshold] <- -1
        chg_dir[is.na(chg_dir)] <- -2
        chg_dir <- array(chg_dir, dim=c(dim(chg_mag)[1], dim(chg_mag)[2], 1))
        return(chg_dir)
    }
    out <- rasterEngine(chg_mag=chg_mag, chg_dir=chg_dir, fun=calc_chg_traj,
                        args=list(chg_threshold=chg_threshold),
                        datatype='INT2S', ...)

    # spatial.tools doesn't properly handle NA values for integer layers, so 
    # they were coded as -2 above - now recode them as NA
    out[out == -2] <- NA

    # spatial.tools can only output the raster package grid format - so output 
    # to a tempfile in that format then copy over to the requested final output 
    # format if a filename was supplied
    if (!missing(filename)) {
        out <- writeRaster(out, filename=filename, overwrite=overwrite, 
                           datatype='INT2S')
    }

    return(out)
}
