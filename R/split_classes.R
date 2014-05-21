#' Split classes in a training dataset using normal mixture modeling
#'
#' This function can be used to aid in classifying spectrally diverse classes 
#' by splitting the input classes into subclasses using a clustering algorithm.  
#' After classification, these subclasses are merged back into their original 
#' parent classes. For example, the training data for an agriculture class 
#' might have both fallow and planted fields in the training data, or fields 
#' planted with different crops that are spectrally dissimilar.  This function 
#' can be used to automatically split the agriculture class into a number of 
#' subclasses.  The classifier is then run on this larger set of classes, and 
#' following classification, these subclasses can all be merged together into a 
#' single overall agriculture class.
#'
#' @export
#' @importFrom mclust Mclust
#' @param train_data a \code{link{pixel_data}} object
#' @param split_levels (optional) a list giving the names of the levels to 
#' split. If missing, all levels will be split.
#' @param verbose whether to report status while running
split_classes <- function(train_data, split_levels, verbose=FALSE) {
    y_reclass <- vector('numeric', nrow(train_data@x))
    if (missing(split_levels)) {
        split_levels <- levels(train_data)
    }
    for (level in split_levels) {
        level_ind <- train_data@y == level
        model <- Mclust(train_data@x[level_ind, ])
        y_reclass[level_ind]  <- paste(train_data@y[level_ind], 
                                       model$classification, sep='_clust')
        if (verbose) print(paste(level, 'split into', model$g, 'classes'))
    }
    y <- factor(y_reclass)
    reclass_mat <- data.frame(split_name=levels(y))
    reclass_mat$split_id <- as.numeric(reclass_mat$split_name)
    reclass_mat$name <- gsub('_clust[0-9]*$', '', reclass_mat$split_name)
    reclass_mat$id <- match(reclass_mat$name, unique(reclass_mat$name))
    return(list(reclass_mat=reclass_mat, y=factor(y_reclass)))
}
