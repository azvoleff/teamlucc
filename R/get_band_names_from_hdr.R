#' Extract the band names listed in an ENVI format header file (.hdr)
#'
#' @export
#' @param hdr_file an ENVI format header file with a .hdr extension
#' @return A /code{list} of band names extracted from the /code{hdr_file}
get_band_names_from_hdr <- function(hdr_file) {
    txt <- readLines(hdr_file)
    line_num <- which(grepl('^band names', txt)) + 1
    band_names <- c()
    for (n in line_num:length(txt)) {
        band_names <- c(band_names, gsub('[,}][[:space:]]*$', '', txt[n]))
        if (grepl('}', txt[n])) {
            break
        }
    }
    return(band_names)
}
