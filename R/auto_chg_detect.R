#' Perform change detection for two Landsat CDR surface reflectance images
#'
#' This image automates the change detection process using the Change Vector 
#' Analysis in Posterior Probability Space (CVAPS) algorithm. The threshold for 
#' change/no-change mapping is determined using Huang's algorithm (see 
#' \code{\link{threshold}}. First the images should be classified using the 
#' \code{auto_classify} function (or any other classification approach that 
#' yields per-pixel probabilities of class membership).
#'
#' @export
#' @param t1_classes cover classes as output from \code{auto_classify_image} 
#' for time 1 image
#' @param t1_probs per class probabilities as output from 
#' \code{auto_classify_image} for time 1 image
#' @param t2_probs per class probabilities as output from 
#' \code{auto_classify_image} for time 2 image
#' @param output_path the path to use for the output
#' @param output_basename the base filename for output files from 
#' \code{auto_chg_detect} (without an extension)
#' @param overwrite whether to overwrite existing files (otherwise an error 
#' will be raised)
#' @param notify notifier to use (defaults to \code{print} function). See the 
#' \code{notifyR} package for one way of sending notifications from R. The 
#' \code{notify} function should accept a string as the only argument.
#' @references Chen, J., X. Chen, X. Cui, and J. Chen. 2011. Change vector 
#' analysis in posterior probability space: a new method for land cover change 
#' detection.  IEEE Geoscience and Remote Sensing Letters 8:317-321.
auto_chg_detect <- function(t1_classes, t1_probs, t2_probs, output_path, 
                            output_basename, overwrite=FALSE, notify=print) {
    if (!file_test("-d", output_path)) {
        stop(paste(output_path, "does not exist"))
    }

    timer <- Track_time(notify)
    timer <- start_timer(timer, label='Change detection')

    ###########################################################################
    # Calculate change magnitude and direction
    ###########################################################################
    timer <- start_timer(timer, label='Change magnitude and direction')

    chg_dir_filename <- file.path(output_path, paste(output_basename, 
                                                     'chgdir.tif', sep='_'))
    chg_dir_image <- chg_dir(t1_probs, t2_probs, filename=chg_dir_filename, 
                             overwrite=overwrite)

    chg_mag_filename <- file.path(output_path, paste(output_basename, 
                                                     'chgmag.tif', sep='_'))
    chg_mag_image <- chg_mag(t1_probs, t2_probs, filename=chg_mag_filename, 
                             overwrite=overwrite)

    timer <- stop_timer(timer, label='Change magnitude and direction')

    ###########################################################################
    # Calculate change trajectories
    ###########################################################################
    timer <- start_timer(timer, label='Change trajectories')
    chg_traj_filename <- file.path(output_path,
                                   paste(output_basename, 'chg_traj.tif', 
                                         sep='_'))

    chg_threshold <- threshold(chg_mag_image, by=.025)

    chg_traj_traj <- chg_traj(t1_classes, chg_mag_image, chg_dir_image, 
                              classnames=levels(t1_classes),
                              chg_threshold=chg_threshold, overwrite=overwrite, 
                              filename=chg_traj_filename)
    timer <- stop_timer(timer, label='Change trajectories')

    timer <- stop_timer(timer, label='Change detection')
}
