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
    cont_table <- table(observed, predicted)
    # Get indices of the diagonal of the cont_table
    diag_indices <- which(diag(nlevels(observed)) == TRUE)
    users_acc <- (cont_table[diag_indices] / colSums(cont_table))  * 100
    prod_acc <- (cont_table[diag_indices] / rowSums(cont_table)) * 100
    overall_acc <- (sum(cont_table[diag_indices]) / sum(cont_table)) * 100
    cont_table <- addmargins(cont_table)
    cont_table <- rbind(cont_table, Users=c(users_acc, NA))
    cont_table <- cbind(cont_table, Producers=c(prod_acc, NA, overall_acc))
    cont_table <- round(cont_table, digits=2)
    dimnames(cont_table) <- list(observed=dimnames(cont_table)[[1]],
                                 predicted=dimnames(cont_table)[[2]])

    return(cont_table)
}
