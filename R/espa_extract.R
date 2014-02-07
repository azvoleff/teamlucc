#' Extract a set of Landsat tarballs into a folder tree organized by image date
#'
#' Each image .tar.gz file will be extracted into a subfolder within 
#' \code{output_folder}. The subfolders will be named according to the year and 
#' Julian date of capture, and the sensor type (LT4, LT5 or LE7 for Landsat 4 
#' TM, Landsat 5 TM, and Landsat 7 ETM+ respectively). For example, 
#' "LT50150531986037-SC20130816144215.tar.gz" would be extracted into a 
#' subfolder named "1986_037_LT5', for 1986, Julian day 37, and Landsat 5 TM.
#'
#' @export
#' @importFrom stringr str_extract
#' @param in_folder Path to a folder of .tar.gz Landsat surface reflectance 
#' images
#' @param out_folder output folder
#' @return used for side effect of unzipping Landsat tarballs
espa_extract <- function(in_folder, out_folder) {
    if (!file_test('-d', in_folder)) {
        stop(paste(in_folder, 'does not exist'))
    }
    if (!file_test('-d', out_folder)) {
        stop(paste(out_folder, 'does not exist'))
    }
    for (zipfile in dir(in_folder, pattern='^.*.tar.gz(ip)?$')) {
        zipfile_path <- file.path(in_folder, zipfile)
        if (!file_test('-f', zipfile_path)) {
            message(paste('Skipping', zipfile, '- not a valid file type.'))
            next
        }
        # Figure out which satellite the image is from
        if (grepl('^LT4', zipfile)) {
            sensor <- 'LT4'
        } else if (grepl('^LT5', zipfile)) {
            sensor <- 'LT5'
        } else if (grepl('^LE7', zipfile)) {
            sensor <- 'LE7'
        } else {
            message(paste('Skipping', zipfile,
                          '- cannot determine sensor from filename.'))
            next
        }
        metadata_string <- str_extract(zipfile, '^((LT4)|(LT5)|(LE7))[0-9]{13}')
        year <- substr(metadata_string, 10, 13)
        julian_day <- substr(metadata_string, 14, 16)
        this_out_folder <- file.path(out_folder,
                                     paste(year, julian_day, sensor, sep='_'))
        if (!file_test('-d', this_out_folder)) {
            dir.create(this_out_folder)
        } else {
            message(paste('Skipping', zipfile, '- output dir', 
                          this_out_folder, 'already exists.'))
            next
        }
        message(paste('***********\nExtracting:\n\t', zipfile, '\nto:\n\t', 
                      this_out_folder))
        ret_code <- untar(zipfile_path, exdir=file.path(this_out_folder))
        if (ret_code != 0) {
            message(paste('WARNING: error extracting', zipfile, '- return code', ret_code))
        }
    }
}
