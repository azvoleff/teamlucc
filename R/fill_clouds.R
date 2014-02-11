# Function to ensure only character variables handed to IDL are quoted
format_IDL_param <- function(varname, varvalue) {
    if (is.character(varvalue)) {
        param <- paste0(varname, '="', varvalue, '"\n')
    } else if (is.list(varvalue)) {
        param <- paste0(varname, '=[')
        for (n in 1:length(varvalue)) {
            if (is.numeric(varvalue[n])) {
                param <- paste0(param, varvalue[n])
            } else {
                param <- paste0(param, '"', varvalue[n], '"')
            }
            if (n != length(varvalue)) {
                param <- paste0(param, ', ')
            }
        }
        param <- paste0(param, ']\n')
    } else {
        param <- paste0(varname, '=', varvalue, '\n')
    }
    return(param)
}

prep_fmask <- function(image_dir) {
    lndsr_regex <- '^lndsr.((LT4)|(LT5)|(LE7)|(LE8))[0-9]{6}[12][0-9]{6}[a-zA-Z]{3}[0-9]{2}'
    lndsr_files <- file.path(image_dir, dir(image_dir, pattern=lndsr_regex))
    fmask_band <- raster(lndsr_files[grepl('fmask_band.envi$', lndsr_files)])

    clouds <- (fmask_band == 2) | (fmask_band == 4)
}

#' Perform heavy cloud filling
#'
#' Calls either of the CLOUD_REMOVE.pro or CLOUD_REMOVE_FAST.pro IDL scripts by 
#' Xiaolin Zhu to fill heavy clouds in a Landsat image.
#'
#' @export
#' @importFrom tools file_path_sans_ext
#' @param img_cloudy the clear image (base image)
#' @param img_clear the cloudy image
#' @param img_cloud_mask cloud mask, with cloud patches assigned unique integer 
#' codes, and all background coded as zero. The \code{ConnCompLabel} function 
#' in the \code{SDMTools} package can automate this process.
#' @param out_name name for output image, or NULL. If null an output name will 
#' be automatically assigned based on \code{img_clear}
#' @param fast if \code{TRUE}, use the CLOUD_REMOVE_FAST.pro script. If 
#' \code{FALSE}, use the CLOUD_REMOVE.pro script.
#' @param num_class set the estimated number of classes in image
#' @param min_pixel the sample size of similar pixels
#' @param extent1 set the range of cloud neighborhood
#' @param DN_min minimum of DN
#' @param DN_max maximum of DN
#' @param patch_long block size when process large image
#' @param idl path to the IDL binary
#' @references Zhu, X., Gao, F., Liu, D., Chen, J., 2012. A modified 
#' neighborhood similar pixel interpolator approach for removing thick clouds 
#' in Landsat images. Geoscience and Remote Sensing Letters, IEEE 9, 521-525.
fill_clouds <- function(img_cloudy, img_clear, img_cloud_mask, out_name=NULL,
                        fast=TRUE, num_class=1, min_pixel=20, extent1=1, 
                        DN_min=0, DN_max=255, patch_long=1000,
                        idl="C:/Program Files/Exelis/IDL83/bin/bin.x86_64/idl.exe") {
    if (fast) {
        script_path <- system.file("idl", "CLOUD_REMOVE_FAST.pro", 
                                   package="teamr")
    } else {
        script_path <- system.file("idl", "CLOUD_REMOVE.pro", 
                                   package="teamr")
    }
    
    if (!(file_test('-x', idl) || file_test('-f', idl))) {
        stop('IDL not found - check "idl" parameter')
    }
    if (!(file_test('-f', img_cloudy))) {
        stop('input file for cloudy image not found- check img_cloudy parameter')
    }
    if (!(file_test('-f', img_clear))) {
        stop('input file for clear image not found- check img_clear parameter')
    }
    if (!(file_test('-f', img_cloud_mask))) {
        stop('input file for cloud mask not found - check img_cloud_mask parameter')
    }

    if (is.null(out_name)) {
        out_name <- paste0(file_path_sans_ext(img_cloudy), '_cloud_remove.envi')
    }

    param_names <- c("cloudy_file", "clear_file", "mask_file", "out_name", 
                     "num_class", "min_pixel", "extent1", "DN_min", "DN_max", 
                     "patch_long")
    param_vals <- list(img_cloudy, img_clear, img_cloud_mask, out_name, 
                       num_class, min_pixel, extent1, DN_min, DN_max, 
                       patch_long)
    idl_params <- mapply(format_IDL_param, param_names, param_vals)
    idl_params <- paste(idl_params, collapse='')

    script_dir <- dirname(script_path)
    idl_script <- tempfile(fileext='.pro')
    idl_cmd <- paste0('CD, "', script_dir, '"\n', idl_params, 'CLOUD_REMOVE,', 
                      paste(param_names, collapse=','), '\nexit')

    f <- file(idl_script, 'wt')
    writeLines(idl_cmd, f)
    close(f)

    idl_out <- system(paste(shQuote(idl), shQuote(idl_script)), intern=TRUE)

    log_file <- paste0(file_path_sans_ext(out_name), '_idllog.txt')
    idl_out <- gsub('\r', '', idl_out)
    f <- file(log_file, 'wt')
    writeLines(idl_out, f) 
    close(f)

    return(brick(out_name))
}
