#' Perform SLC-off gap fill of Landsat 7 ETM+ image
#'
#' Calls GNSPI.pro IDL script by Xiaolin Zhu to fill gaps in SLC-off Landsat 7 
#' ETM+ image.
#'
#' @importFrom tools file_path_sans_ext
#' @param slc_off_file the SLC-off Landsat 7 file to gap fill
#' @param input_file the first file to use to fill in the gaps
#' @param timeseries_files a timeseries of files to use as additional inputs to 
#' the gap fill algorithm
#' @param out_base path and base filename for the output file. The script will 
#' save the output files by appending "_GNSPI.envi" and 
#' "_GNSPI_uncertainty.envi" to this base filename.
#' @param sample_size the sample size of sample pixels
#' @param size_wind the maximum window size
#' @param class_num the estimated number of classes
#' @param DN_min the minimum DN value of the image
#' @param DN_max the maximum DN value of the image
#' @param patch_long the size of block, to process whole ETM scene, set to 1000
#' @param idl path to the IDL binary
#' @export
#' @references Zhu, X., Liu, D., Chen, J., 2012. A new geostatistical approach 
#' for filling gaps in Landsat ETM+ SLC-off images. Remote Sensing of 
#' Environment 124, 49-60.
fill_gaps <- function(slc_off_file, input_file, timeseries_files, 
                      out_base=NULL, sample_size=20, size_wind=12, class_num=4, 
                      DN_min=0.0, DN_max=1.0, patch_long=1000,
                      idl="C:/Program Files/Exelis/IDL83/bin/bin.x86_64/idl.exe") {
    script_path <- system.file("idl", "GNSPI.pro", package="teamr")

    if (!(file_test('-x', idl) || file_test('-f', idl))) {
        stop('IDL not found - check "idl" parameter')
    }
    if (!(file_test('-f', slc_off_file))) {
        stop('SLC-off file not found- check slc_off_file parameter')
    }
    if (!(file_test('-f', input_file))) {
        stop('input file for gap fill not found- check input_file parameter')
    }
    for (n in 1:length(timeseries_files)) {
        if (!(file_test('-f', timeseries_files[[n]]))) {
            stop(paste0('timeseries file "', timeseries_files[[n]], '" not found'))
        }
    }

    if (is.null(out_base)) {
        out_base <- paste0(file_path_sans_ext(slc_off_file))
    }

    temp_dir <- tempdir()

    param_vals <- list(slc_off_file, input_file, timeseries_files,
                       out_base, sample_size, size_wind, class_num,
                       DN_min, DN_max, patch_long, temp_dir)
    param_names <- list('slc_off_file', 'input_file', 'timeseries_files', 
                        'out_base', 'sample_size', 'size_wind', 'class_num', 
                        'DN_min', 'DN_max', 'patch_long', 'temp_dir')
    idl_params <- mapply(format_IDL_param, param_names, param_vals)
    idl_params <- paste(idl_params, collapse='')

    script_dir <- dirname(script_path)
    idl_script <- tempfile(fileext='.pro')
    idl_cmd <- paste0('CD, "', script_dir, '"\n', idl_params, 'GNSPI,', 
                      paste(param_names, collapse=','), '\nexit')

    f <- file(idl_script, 'wt')
    writeLines(idl_cmd, f)
    close(f)

    idl_out <- system(paste(shQuote(idl), shQuote(idl_script)), intern=TRUE)

    log_file <- paste0(out_base, '_GNSPI_idllog.txt')
    idl_out <- gsub('\r', '', idl_out)
    f <- file(log_file, 'wt')
    writeLines(idl_out, f) 
    close(f)

    return(brick(paste0(out_base, '_GNSPI.envi')))
}
