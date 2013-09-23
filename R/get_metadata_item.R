#' Extract a metadata item from a metadata file in team format 
#'
#' @export
#' @param x a metadata file in team format (ending in .txt)
#' @param item a string giving the name of the metadata item to extract
#' @return The metadata item (as a string)
#' @examples
#' # Load an example metadata file
#' L5TSR_1986_file <- system.file('extdata/L5TSR_1986.dat', package='teamr')
#' x <- extension(L5TSR_1986_file, 'txt')
#'
#' # Extract the sun angle elevation from the metadata file. Note that
#' # get_metadata_item returns a string, so the item must be converted using
#' # as.numeric
#' sunelev <- as.numeric(get_metadata_item(x, 'SolarZenith'))
get_metadata_item <- function(x, item) {
    if (!file.exists(x)) {
        stop(paste('Could not find metadata file', x))
    }
    metadata <- read.csv(x, stringsAsFactors=FALSE)
    rownum <- which(grepl(item, metadata$item))
    return(metadata$value[rownum])
}
