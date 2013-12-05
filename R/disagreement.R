#' Calculate quantity and allocation disagreement for two maps
#'
#' This function calculate quantity disagreement \code{Q} and allocation 
#' disagreement \code{A} for use in accuracy assessment and map comparison.  
#' These quantities are calculated based on Pontius and Millones (2011).
#'
#' @export
#' @param cont_table an estimated population matrix (contingency table), or, if 
#' \code{ref_map} is provided, an observed sample matrix as calculated by 
#' \code{accuracy}.
#' @param ref_map (optional) a baseline reference map as a \code{RasterLayer} 
#' to be used for converting the observed sample matrix (\code{cont_table}) to 
#' an estimated population matrix (using equation 1 in Pontius and Millones 
#' 2011).
#' @param comp_map (optional) a comparison map as a \code{RasterLayer}, to be 
#' used for converting the observed sample matrix (\code{cont_table}) to an 
#' estimated population matrix (using equation 1 in Pontius and Millones 2011).  
#' If ref_map is provided to \code{disagreement}, but comp_map is not provided, 
#' then \code{disagreement} will assume that the class proportions in ref_map 
#' and comp_map are equal.  Testing the effects of this assumption is highly 
#' recommended.
#' @return list with two entries: Q (quantity disagreement) and A (allocation 
#' disagreement)
#' @references Pontius, R. G., and M. Millones. 2011. Death to Kappa: birth of 
#' quantity disagreement and allocation disagreement for accuracy assessment.  
#' International Journal of Remote Sensing 32:4407-4429.
disagreement <- function(cont_t, ref_map=NULL, comp_map=NULL) {
    if (is.null(ref_map) & (!is.null(comp_map))) {
        stop('if comp_map is provided, ref_map must also be provided')
    } else if (!is.null(ref_map)) {
        # cont_table is an observed sample matrix - convert it to an estimated 
        # population matrix
    }
}
