#' Extract a metadata item from a metadata file in GDAL PAM format
#'
#' GDAL PAM format metadata files end in ".aux.xml".
#'
#' @export
#' @importFrom raster extension
#' @importFrom XML xmlInternalTreeParse xpathApply xmlValue
#' @param x an image file that has an accompanying GDAL PAM format metadata 
#' file (ending in .aux.xml)
#' @param key a string giving the name of the metadata item to extract
#' @return The metadata item (as a string)
get_metadata_item <- function(x, key) {
    metadata_file <- paste0(x, '.aux.xml')
    if (!file.exists(metadata_file)) {
        stop(paste('Could not find metadata file', metadata_file))
    }
    doc <- xmlInternalTreeParse(metadata_file)
    xpath_exp <- paste0("//MDI[@key='", key, "']")
    value <- unlist(xpathApply(doc, xpath_exp, xmlValue))
    if (length(value) > 1) {
        stop('multiple elements found')
    }
    if (length(value) == 0) {
        stop('no elements found')
    }
    return(value)
}
