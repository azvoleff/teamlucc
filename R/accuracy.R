#' Calculate statistics on classification accuracy
#'
#' Used for calculating a contingency table with user's, producer's, and 
#' overall accuracies for an image classification.
#'
#' @export
#' @param model a classification model with a \code{predict} method
#' @param test_data a training/testing dataset as output by the 
#' \code{extract_training_data} function. If a 'Training' column is included 
#' \code{accuracy} will use only the data not used in training for evaluating 
#' model accuracy.
#' @return list of accuracy statistics
#' @examples
accuracy <- function(model, test_data) {
    if (!('Training' %in% names(test_data))) {
        warning('no Training column found - assuming none of "test_data" was used for model training')
    } else {
        if (sum(test_data$Training) == 0) {
            stop('cannot conduct accuracy assessment without independent testing data')
        }
        test_data <- test_data[!test_data$Training, ]
    }
    observed <- test_data$y
    predicted <- predict(model, test_data)
    ct <- table(observed, predicted)
    # Get indices of the diagonal of the ct
    diag_indices <- which(diag(nlevels(observed)) == TRUE)
    users_acc <- (ct[diag_indices] / colSums(ct))  * 100
    prod_acc <- (ct[diag_indices] / rowSums(ct)) * 100
    overall_acc <- (sum(ct[diag_indices]) / sum(ct)) * 100
    ct <- addmargins(ct)
    ct <- rbind(ct, Users=c(users_acc, NA))
    ct <- cbind(ct, Producers=c(prod_acc, NA, overall_acc))
    ct <- round(ct, digits=2)
    dimnames(ct) <- list(observed=dimnames(ct)[[1]],
                                 predicted=dimnames(ct)[[2]])
    return(ct)
}
