#' Calculate a contingency table using the composite operator
#'
#' This function calculates a cross tabulation for a map comparison using the 
#' composite operator recommended by Pontius and Cheuk (2006).
#'
#' @export
#' @param model a classification model with a \code{predict} method
#' @param test_data a training/testing dataset as output by the 
#' \code{extract_training_data} function. If a 'Training' column is included 
#' \code{accuracy} will use only the data not used in training for evaluating 
#' model accuracy.
#' @return contingency table
#' @references Pontius, R. G., and M. Millones. 2011. Death to Kappa: birth of 
#' quantity disagreement and allocation disagreement for accuracy assessment.  
#' International Journal of Remote Sensing 32:4407-4429.
compcont <- function(cont_t, ref_map=NULL, comp_map=NULL) {
    stop('compcont is not yet finished')
    #TODO: finish coding this function
}
