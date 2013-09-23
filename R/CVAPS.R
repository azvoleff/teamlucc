#' Change Vector Analysis in Posterior Probability Space
#'
#' @export
#' @param svm_t0 Support Vector Machine for time 0
#' @param x_t0 Raster object with predictor layers for time 0
#' @param svm_t1 Support Vector Machine for time 1
#' @param x_t1 Raster object with predictor layers for time 1
#' @references Chen, J., X. Chen, X. Cui, and J. Chen. 2011. Change vector 
#' analysis in posterior probability space: a new method for land cover change 
#' detection. IEEE Geoscience and Remote Sensing Letters 8:317–321.
CVAPS <- function(svm_t0, x_t0, svm_t1, x_t1) {
    # Check that x_t0 and x_t1 are matched images (identical extents and 
    # identical projections)
    
    # Check that svm_t0 and svm_t1 have the same classes
    
    # Get number of classes from coefs matrix
    n_classes <- svm_t0@coefs

    pred_t0 <- predict(best.tune(svm_t0), newdata=x_t0, probability=TRUE)
    pred_t1 <- predict(best.tune(svm_t1), newdata=x_t1, probability=TRUE)
    t0_prob <- pred_t0@probabilities
    t1_prob <- pred_t1@probabilities

    # May need to run getValues first before running the below calculations

    # Calculate change magnitude (eqn 3 in Chen 2011)
    dP <- t1_prob - t0_prob
    chg_mag <- sqrt(rowSums((t1_prob - t0_prob)^2))

    # Calculate change vector (eqn 3 in Chen 2011). Use dot product of the base 
    # change vectors Ea,b,. Each base change vector represents complete change 
    # from class a to class b.
    Eab <- diag(n_classes)
    chg_type <- apply(dP %*% Eab, 1, function(r) which(r == max(r)))

    # Now do double-window flexible pace search
    
}
