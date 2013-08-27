get_metadata_item <- function(metadatafile, item) {
    if (!file.exists(metadatafile)) {
        stop(paste('Could not find metadata file', metadatafile))
    }
    metadata <- read.csv(metadatafile, stringsAsFactors=FALSE)
    rownum <- which(grepl(item, metadata$item))
    return(metadata$value[rownum])
}
