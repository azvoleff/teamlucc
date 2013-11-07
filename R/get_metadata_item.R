#' Extract a metadata item from a metadata file in team format 
#'
#' @export
#' @param x a metadata file in team format (ending in .txt)
#' @param item a string giving the name of the metadata item to extract
#' @return The metadata item (as a string)
get_metadata_item <- function(x, item) {
    if (!file.exists(x)) {
        stop(paste('Could not find metadata file', x))
    }
    metadata <- read.csv(x, stringsAsFactors=FALSE)
    rownum <- which(grepl(item, metadata$item))
    return(metadata$value[rownum])
}
