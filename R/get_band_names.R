get_band_names <- function(hdr_file) {
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

for (n in 1:nrow(image_list)) {
    in_prefix <- sub('.hdf$', '', image_list[n, ]$file_path)
    out_prefix <- sub('orig', 'proc', in_prefix)
    print(paste('Processing', in_prefix))

    fc <- brick(paste(out_prefix, '.bsq', sep=''))
    names(fc) <- get_band_names(paste(out_prefix, '.hdr', sep=''))

    fc <- stack(paste(in_file_prefix, '_band4.bsq', sep=''),
                         paste(in_file_prefix, '_band3.bsq', sep=''),
                         paste(in_file_prefix, '_band2.bsq', sep=''))

}
