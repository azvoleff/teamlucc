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
team_chg_detect <- function(t1_predictions, t2_predictions) {
    ################################################################################
    # Calculate change magnitude and direction
    ################################################################################
    print('Calculating change magnitude/direction...')
    trackTime(action='start')
    t0_probs <- brick('L5TSR_1986_probs.grd')
    t1_probs <- brick('L5TSR_2001_probs.grd')
    chg_dir_image <- chg_dir(t0_probs, t1_probs, filename='L5TSR_1986_to_2001_chgdir')
    chg_mag_image <- chg_mag(t0_probs, t1_probs, filename='L5TSR_1986_to_2001_chgmag')
    trackTime()

    ################################################################################
    # Calculate change trajectories
    ################################################################################

    sfQuickStop()
    trackTime('total_time')
}
