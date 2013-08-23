library(rgdal)
library(raster)
library(landsat)

PLOT_WIDTH <- 3.25
PLOT_HEIGHT <- 3.25
DPI <- 300

image_list <- read.csv('H:/Data/TEAM/VB/Rasters/Landsat/image_list_test.csv')
out_folder <- 'H:/Data/TEAM/VB/Rasters/Landsat'

get_band_names <- function(hdr_file) {
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
