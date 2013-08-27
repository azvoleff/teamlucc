get_band_names_from_hdr <- function(hdr_file) {
    require(raster)
    txt <- readLines(hdr_file)
    line_num <- which(grepl('^band names', a)) + 1
    band_names <- c()
    for (n in line_num:length(txt)) {
        band_names <- c(band_names, gsub('[,}][[:space:]]*$', '', txt[n]))
        if (grepl('}', txt[n])) {
            break
        }
    }
    return(band_names)
}
