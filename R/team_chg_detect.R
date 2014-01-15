#' Perform change detection for two Landsat CDR surface reflectance images
#'
#' First the images should be classified using the \code{team_classify} 
#' function.
#'
#' @export
#' @param t1_predictions output from \code{team_classify_image} for time 1 
#' image
#' @param t2_predictions output from \code{team_classify_image} for time 2 
#' image
#' @param notify notifier to use (defaults to \code{print} function). See the 
#' \code{notifyR} package for one way of sending notifications from R. The 
#' \code{notify} function should accept a string as the only argument.
team_chg_detect <- function(t1_predictions, t2_predictions, notify) {
    timer <- Track_time(notify)

    timer <- start_timer(timer, label=paste('Change detection', image_basename))

    ################################################################################
    # Calculate change magnitude and direction
    ################################################################################
    timer <- start_timer(timer, label=paste('Change magnitude and direction', image_basename))
    t0_probs <- brick('L5TSR_1986_probs.grd')
    t1_probs <- brick('L5TSR_2001_probs.grd')
    chg_dir_image <- chg_dir(t0_probs, t1_probs, filename='L5TSR_1986_to_2001_chgdir')
    chg_mag_image <- chg_mag(t0_probs, t1_probs, filename='L5TSR_1986_to_2001_chgmag')
    timer <- stop_timer(timer, label=paste('Change magnitude and direction', image_basename))

    ################################################################################
    # Calculate change trajectories
    ################################################################################
    timer <- start_timer(timer, label=paste('Change trajectories', image_basename))
    timer <- stop_timer(timer, label=paste('Change trajectories', image_basename))

    sfQuickStop()

    timer <- stop_timer(timer, label=paste('Change detection', image_basename))
}
