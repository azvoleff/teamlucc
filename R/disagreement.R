#' Calculate quantity and allocation disagreement for two maps
#'
#' This function calculate quantity disagreement \code{Q} and allocation 
#' disagreement \code{A} for use in accuracy assessment and map comparison.  
#' These quantities are calculated based on Pontius and Millones (2011).
#'
#' @export
#' @param an observed sample matrix as calculated by \code{accuracy}
#' @param pop (optional) either 1) a list of length equal to the number of 
#' classes in ct giving the total number of pixels in the entire
#' population for each class, or 2) the predicted cover map as a 
#' \code{RasterLayer}, from which the population frequencies will be derived.
#' @return list with two entries: Q (quantity disagreement) and A (allocation 
#' disagreement)
#' @references Pontius, R. G., and M. Millones. 2011. Death to Kappa: birth of 
#' quantity disagreement and allocation disagreement for accuracy assessment.  
#' International Journal of Remote Sensing 32:4407-4429.
disagreement <- function(ct, pop=NULL) {
    if (nrow(ct) != ncol(ct)) {
        stop('ct must be square')
    }
    if ('RasterLayer' %in% class(pop)) {
        pop <- freq(pop)[1]
        if (length(pop) != nrow(ct)) {
            stop('number of classes in pop must be equal to nrow(ct)')
        }
    } else if (is.null(pop)) {
        warning('pop was not provided - assuming sample frequencies equal population frequencies')
        pop <- rowSums(ct)
    } else if (class(pop) %in% c('integer', 'numeric')) {
        if (length(pop) != nrow(ct)) {
            stop('length(pop) must be equal to nrow(ct)')
        }
    } else { 
        stop('pop must be a numeric vector, integer vector, RasterLayer, or NULL')
    }
    # Below uses the notation of Pontius and Millones (2011)
    nijsum <- matrix(rowSums(ct), nrow=nrow(ct), ncol=ncol(ct))
    Ni <- matrix(pop, nrow=nrow(ct), ncol=ncol(ct))
    pop_ct <- (ct / nijsum) * (Ni / sum(pop))

    # Calculate quantity disagreement (Pontius and Millones, 2011, eqns 2-3)
    qg_mat = abs(rowSums(pop_ct) - colSums(pop_ct))
    Q = sum(qg_mat) / 2

    # Calculate allocation disagreement (Pontius and Millones, 2011, eqns 4-5)
    diag_indices <- which(diag(nrow(pop_ct)) == TRUE)
    ag_mat = 2 * apply(cbind(rowSums(pop_ct) - pop_ct[diag_indices],
                             colSums(pop_ct) - pop_ct[diag_indices]), 1, min)
    A = sum(ag_mat) / 2

    return(list(Q=Q, A=A))
}
