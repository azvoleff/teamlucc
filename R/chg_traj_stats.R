#' Calculate change-trajectory statistics
#'
#' @export
#' @param traj a list (as output by \code{chg_traj} with two elements: traj_lut 
#' (a lookup table of change trajectory codes) and chg_traj (a 
#' \code{RasterLayer} of change trajectory codes.
#' @return a \code{data.frame} object with change trajectory statistics
#' @examples
#' #TODO: Add examples
chg_traj_stats <- function(traj) {
    chg_table <- table(getValues(traj$chg_traj))
    summ_table <- data.frame(Traj_Code=traj$traj_lut$Code,
                             Trajectory=paste(traj$traj_lut$t0_name, traj$traj_lut$t1_name, 
                                              sep='-'))
    summ_table <- cbind(summ_table, n_pixels=chg_table[match(row.names(chg_table), summ_table$Traj_Code)])
    row.names(summ_table) <- NULL
    summ_table$Frac_Chg <- summ_table$n_pixels / sum(summ_table$n_pixels)
    summ_table$Frac_Tot <- summ_table$n_pixels / length(traj$chg_traj)
    return(summ_table)
}
