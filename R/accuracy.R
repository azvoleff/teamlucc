setClass("accuracy",
         representation(genome="table", # a single string
                        snpid="character", # a character vector of length N
                        chrom="character", # a character vector of length N
                        pos="integer" # an integer vector of length N
                        )
         )

.calc_pop_ct <- function(ct, pop) {
    # Below uses the notation of Pontius and Millones (2011)
    nijsum <- matrix(rowSums(ct), nrow=nrow(ct), ncol=ncol(ct))
    Ni <- matrix(pop, nrow=nrow(ct), ncol=ncol(ct))
    # pop_ct is the population contigency table
    return((ct / nijsum) * (Ni / sum(pop)))
}

.calc_Q <- function(pop_ct) {
    # Calculate quantity disagreement (Pontius and Millones, 2011, eqns 2-3)
    qg_mat = abs(rowSums(pop_ct) - colSums(pop_ct))
    return(sum(qg_mat) / 2)
}

.calc_A <- function(pop_ct) {
    # Calculate allocation disagreement (Pontius and Millones, 2011, eqns 4-5)
    diag_indices <- which(diag(nrow(pop_ct)) == TRUE)
    ag_mat = 2 * apply(cbind(rowSums(pop_ct) - pop_ct[diag_indices],
                             colSums(pop_ct) - pop_ct[diag_indices]), 1, min)
    return(sum(ag_mat) / 2)
}

.add_ct_margins <- function(ct) {
    # Get column-major indices of the diagonal of the ct
    diag_indices <- which(diag(nrow(ct)) == TRUE)
    users_acc <- ct[diag_indices] / colSums(ct)
    prod_acc <- ct[diag_indices] / rowSums(ct)
    overall_acc <- sum(ct[diag_indices]) / sum(ct)
    ct <- addmargins(ct)
    ct <- rbind(ct, Producers=c(users_acc, NA))
    ct <- cbind(ct, Users=c(prod_acc, NA, overall_acc))
    ct <- round(ct, digits=4)
    dimnames(ct) <- list(predicted=dimnames(ct)[[1]],
                         observed=dimnames(ct)[[2]])
    return(ct)
}

#' Calculate statistics summarizing classification accuracy
#'
#' Calculates a contingency table and various statistics for use in image 
#' classification accuracy assessment and map comparison. Contingency table 
#' includes user's, producer's, and overall accuracies for an image 
#' classification, and quantity disagreement \code{Q} and allocation 
#' disagreement \code{A}. Q and A are calculated based on Pontius and Millones 
#' (2011). 95% confidence intervals for the user's, producer's, and overall 
#' accuracies are calculated as in Olofsson et al. (2013).
#'
#' To avoid bias due to the use of a sample contingency table, the contingency 
#' table can be converted to a population contingency table, if the variable 
#' 'pop' is provided. For an accuracy assessment based on a random sample, 
#' 'pop' does not need to be provided.
#'
#' @export
#' @param model a classification model with a \code{predict} method
#' @param test_data a training/testing dataset as output by the 
#' \code{extract_training_data} function. If a 'Training' column is included, 
#' \code{accuracy} will use only the data not used in training for evaluating 
#' model accuracy. If test_data is NULL, \code{accuracy} will try to use the 
#' trainingData included in the model object (this will only work if the model 
#' object is of class \code{train} from the \code{caret} package.
#' @param pop (optional) used to convert from sample matrix to population 
#' matrix as in Pontius and Millones 2011.Can be: 1) NULL, in which case the 
#' sample frequencies will be used as estimates of the population frequencies, 
#' 2) a list of length equal to the number of classes in the map giving the 
#' total number of pixels in the entire population for each class, or 3) the 
#' predicted cover map from \code{model} as a \code{RasterLayer}, from which 
#' the population frequencies will be tabulated.
#' @return list with four elements: ct (sample contingency table), pop_ct 
#' (population contingency table), Q (quantity disagreement), and A (allocation 
#' disagreement).
#' @references Pontius, R. G., and M. Millones. 2011. Death to Kappa: birth of 
#' quantity disagreement and allocation disagreement for accuracy assessment.  
#' International Journal of Remote Sensing 32:4407-4429.
#' Olofsson, P., G. M. Foody, S. V. Stehman, and C. E. Woodcock.  2013. Making 
#' better use of accuracy data in land change studies: Estimating accuracy and 
#' area and quantifying uncertainty using stratified estimation.  Remote 
#' Sensing of Environment 129:122-131.
#' @examples
#' accuracy(classified_LT5SR_1986$model)
accuracy <- function(model, test_data=NULL, pop=NULL) {
    if (('train' %in% class(model)) && is.null(test_data)) {
        test_data <- model$trainingData
        names(test_data)[names(test_data) == '.outcome'] <- 'y'
    }
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
    # ct is the sample contigency table
    ct <- table(predicted, observed)

    if ('RasterLayer' %in% class(pop)) {
        pop <- freq(pop, useNA='no')[, 2]
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

    pop_ct <- .calc_pop_ct(ct, pop)
    Q <- .calc_Q(pop_ct)
    A <- .calc_A(pop_ct)

    return(list(ct=.add_ct_margins(ct), pop_ct=.add_ct_margins(pop_ct),
                Q=Q, A=A))
}
