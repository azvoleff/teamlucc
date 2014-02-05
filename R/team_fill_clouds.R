# Function to ensure only character variables handed to IDL are quoted
format_IDL_arg <- function(varname, varvalue) {
    if (is.character(varvalue)) {
        paste0(varname, '="', varvalue, '"\n')
    } else {
        paste0(varname, '=', varvalue, '\n')
    }
}

#' Perform heavy cloud filling
#'
#' Calls CLOUD_REMOVE.pro IDL script by Zhu Xiaolin to fill heavy clouds in 
#' Landsat image.
#'
#' @export
#' @importFrom tools file_path_sans_ext
#' @param img_cloudy the clear image (base image)
#' @param img_clear the cloudy image
#' @param img_cloud_mask cloud mask, with clouded areas set to
#' @param script_path the path to the CLOUD_REMOVE.pro file
#' @param num_class set the estimated number of classes in image
#' @param extent1 set the range of cloud neighborhood
#' @param DN_min minimum of DN
#' @param DN_max maximum of DN
#' @param patch_long block size when process large image
#' @param idl the path to the IDL binary
#' @references Zhu, X., Gao, F., Liu, D., Chen, J., 2012. A modified 
#' neighborhood similar pixel interpolator approach for removing thick clouds 
#' in Landsat images. Geoscience and Remote Sensing Letters, IEEE 9, 521-525.
team_fill_clouds <- function(img_cloudy, img_clear, img_cloud_mask, out_name=NULL,
                             script_path='CLOUD_REMOVE.sav', num_class=1, 
                             extent1=1, DN_min=0, DN_max=255, patch_long=1000,
                             idl="C:/Program Files/Exelis/IDL83/bin/bin.x86_64/idl.exe") {
    # img_cloudy <- 'C:/Users/azvoleff/Code/TEAM/team_IDL_code/cloud_removal/test data/L20080724_cloudy'
    # img_clear <- 'C:/Users/azvoleff/Code/TEAM/team_IDL_code/cloud_removal/test data/L20080606'
    # img_cloud_mask <- 'C:/Users/azvoleff/Code/TEAM/team_IDL_code/cloud_removal/test data/cloud_mask'
    # script_path <- 'C:/Users/azvoleff/Code/TEAM/team_IDL_code/cloud_removal/CLOUD_REMOVE.pro'
    # out_name=NULL
    # num_class <- 1
    # extent1 <- 1
    # DN_min <- 0
    # DN_max <- 255
    # patch_long <- 1000
    # idl <- "C:/Program Files/Exelis/IDL83/bin/bin.x86_64/idl.exe"

    if (is.null(out_name)) {
        out_name <- paste0(file_path_sans_ext(img_cloudy), '_cloud_remove.envi')
    }

    varnames <- c("cloudy_file", "clear_file", "mask_file", "out_name", 
                  "num_class", "extent1", "DN_min", "DN_max", "patch_long")
    varvalues <- list(img_cloudy, img_clear, img_cloud_mask, out_name, 
                      num_class, extent1, DN_min, DN_max, patch_long)

    script_vars <- mapply(format_IDL_arg, varnames, varvalues)
    script_vars <- paste(script_vars, collapse='')

    script_dir <- dirname(script_path)
    idl_script <- tempfile(fileext='.pro')
    idl_cmd <- paste0('CD, "', script_dir, '"\n', script_vars, 'CLOUD_REMOVE,', paste(varnames, collapse=','), '\nexit')

    f <- file(idl_script, 'wt')
    writeLines(idl_cmd, f)
    close(f)

    idl_out <- system(paste(shQuote(idl), shQuote(idl_script)), intern=TRUE)

    log_file <- paste0(file_path_sans_ext(out_name), '_log.txt')
    idl_out <- gsub('\r', '', idl_out)
    f <- file(log_file, 'wt')
    writeLines(idl_out, f) 
    close(f)

    return(brick(out_name))
}
