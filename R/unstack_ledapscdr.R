#' Convert Landsat CDR images from HDF4 to ENVI format
#'
#' This function converts a Landsat surface reflectance (SR) image from the 
#' Landsat Climate Data Record (CDR) archive into a series of single band 
#' images in ENVI format.
#'
#' Note that this function requires a local GDAL installation. See 
#' http://trac.osgeo.org/gdal/wiki/DownloadingGdalBinaries or 
#' http://trac.osgeo.org/osgeo4w/ to download the appropriate installer for 
#' your operating system.
#'
#' @export
#' @importFrom tools file_path_sans_ext
#' @param x input HDF4 file
#' @param out_folder output folder (if \code{NULL}, defaults to input folder)
#' @param overwrite whether to overwrite existing files
#' @param rmhdf whether to remove hdf files after unstacking them
#' @return nothing (used for side effect of converting Landsat CDR HDF files)
#' @examples
#' \dontrun{
#' # Unstack files downloaded from CDR:
#' unstack_ledapscdr('lndsr.LT50150531986021XXX03.hdf')
#' 
#' # Unstack files downloaded from CDR, overwriting any existing files, and 
#' # deleting original HDF files after unstacking:
#' unstack_ledapscdr('lndsr.LT50150531986021XXX03.hdf', overwrite=TRUE, 
#'                   rmhdf=TRUE)
#' }
unstack_ledapscdr <- function(x, output_folder=NULL, overwrite=FALSE, 
                              rmhdf=FALSE) {
    if (!require(gdalUtils)) {
        stop('unstack_ledapscdr requires the "gdalUtils" package')
    }

    if (is.null(output_folder)) {
        output_folder <- dirname(x)
    }

    if ((!file_test('-f', x)) | (tolower(extension(x)) != '.hdf')) {
        stop('x must be an existing file ending in ".hdf"')
    }

    if (!file_test('-d', output_folder)) {
        stop('output_folder must be a directory')
    }
    out_basename <- file_path_sans_ext(basename(x))

    sds <- get_subdatasets(x)
    loc <- regexpr('[a-zA-Z0-9_-]*$', sds)
    for (n in 1:length(sds)) {
        start_char <- loc[n]
        stop_char <- start_char + attr(loc, 'match.length')[n]
        band_name <- substr(sds[[n]], start_char, stop_char)
        this_out <- paste0(file.path(output_folder, out_basename), '_', band_name, '.envi')
        if (file.exists(this_out)) {
            if (overwrite) {
                unlink(this_out)
                if (file.exists(extension(this_out, 'hdr'))) 
                    unlink(extension(this_out, 'hdr'))
                if (file.exists(extension(this_out, 'envi.aux.xml'))) 
                    unlink(extension(this_out, 'envi.aux.xml'))
                if (file.exists(extension(this_out, 'envi.enp'))) 
                    unlink(extension(this_out, 'envi.enp'))
            } else {
                warning(paste(this_out, 'already exists - skipping file'))
                next
            }
        }
        out_rast <- gdal_translate(x, of="ENVI", sd_index=n, this_out, 
                                   outRaster=TRUE)
    }
    if (rmhdf) {
        if (file.exists(extension(x, 'hdf'))) 
            unlink(extension(x, 'hdf'))
        if (file.exists(extension(x, 'hdr'))) 
            unlink(extension(x, 'hdr'))
        if (file.exists(paste0(file_path_sans_ext(x), '.hdf.hdr')))
            unlink(paste0(file_path_sans_ext(x), '.hdf.hdr'))
        if (file.exists(extension(x, 'txt'))) 
            unlink(extension(x, 'txt'))
    }
}
