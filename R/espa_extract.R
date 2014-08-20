#' Extract a set of Landsat tarballs into a folder tree organized by image date
#'
#' Each image .tar.gz file will be extracted into a subfolder within 
#' \code{output_folder}. The subfolders will be named according to the year and 
#' Julian date of capture, and the sensor type (LT4, LT5 or LE7 for Landsat 4 
#' TM, Landsat 5 TM, and Landsat 7 ETM+ respectively). For example, 
#' "LT50150531986037-SC20130816144215.tar.gz" would be extracted into a 
#' subfolder named "1986_037_LT5', for 1986, Julian day 37, and Landsat 5 TM.
#'
#' Zip files for images from sensors not included in \code{sensors}, from 
#' path/rows not included in \code{pathrows} or with acquisition dates
#' outside of the period defined by \code{start_date} and \code{end_date} will 
#' be ignored.
#'
#' @export
#' @importFrom stringr str_extract
#' @param in_folder Path to a folder of .tar.gz Landsat surface reflectance 
#' images
#' @param out_folder output folder
#' @param pathrows a list of paths/rows to include. Each path/row should be 
#' specified as a six digit string. For example, path 231, row 62 would be 
#' specified as "231062".
#' @param start_date start date of period from which images will be extracted
#' to (as \code{Date} object).
#' @param end_date end date of period from which images will be extracted
#' to (as \code{Date} object)
#' @param sensors a list of the sensors to include (can be any of "LT4", "LT5", 
#' "LE7", or "LC8")
#' @return nothing (used for side effect of unzipping Landsat CDR tarballs)
#' @examples
#' \dontrun{
#' # Don't filter:
#' espa_extract('D:/Landsat_Originals', 'D:/Landsat_Out')
#'
#' # Filter by start and end date:
#' start_date <- as.Date('2010/1/1')
#' end_date <- as.Date('2010/12/31')
#' espa_extract('D:/Landsat_Originals', 'D:/Landsat_Out',
#'              start_date=start_date, end_date=end_date)
#'
#' # Filter by start and end date, sensor, and pathrow:
#' espa_extract('D:/Landsat_Originals', 'D:/Landsat_Out', 
#'              start_date=start_date, end_date=end_date, sensors='LE7',
#'              pathrows='231062')
#' }
espa_extract <- function(in_folder, out_folder, pathrows=NULL, start_date=NULL, 
                         end_date=NULL, sensors=NULL) {
    if (!file_test('-d', in_folder)) {
        stop(paste(in_folder, 'does not exist'))
    }
    if (!file_test('-d', out_folder)) {
        stop(paste(out_folder, 'does not exist'))
    }

    zipfiles <- dir(in_folder, pattern='^.*.tar.gz(ip)?$')

    # Filter by date
    img_dates <- as.Date(gsub('-', '', str_extract(zipfiles, '[0-9]{7}-')), '%Y%j')
    if (!is.null(start_date)) {
        stopifnot(class(start_date) == 'Date')
        inc_dates <- which(img_dates >= start_date)
        zipfiles <- zipfiles[inc_dates]
        img_dates <- img_dates[inc_dates]
    }
    if (!is.null(end_date)) {
        stopifnot(class(end_date) == 'Date')
        inc_dates <- which(img_dates < end_date)
        zipfiles <- zipfiles[inc_dates]
        img_dates <- img_dates[inc_dates]
    }

    # Filter by pathrow
    img_pathrows <- gsub('(LT[45])|(LE7)|(LC8)', '', str_extract(zipfiles, '((LT[45])|(LE7)|(LC8))[0-9]{6}'))
    if (!is.null(pathrows)) {
        stopifnot(!is.na(str_extract(pathrows, '[0-9]{6}')))
        inc_pathrows <- img_pathrows %in% pathrows
        zipfiles <- zipfiles[inc_pathrows]
        img_pathrows <- img_pathrows[inc_pathrows]
        img_dates <- img_dates[inc_pathrows]
    }

    # Filter by sensor
    img_sensors <- str_extract(zipfiles, '^((LT[45])|(LE7)|(LC8))')
    if (!is.null(sensors)) {
        stopifnot(!is.na(str_extract(sensors, '^((LT[45])|(LE7)|(LC8))$')))
        inc_sensors <- img_sensors %in% sensors
        zipfiles <- zipfiles[inc_sensors]
        img_pathrows <- img_pathrows[inc_sensors]
        img_sensors <- img_sensors[inc_sensors]
        img_dates <- img_dates[inc_sensors]
    }

    img_paths <- str_extract(img_pathrows, '^[0-9]{3}')
    img_rows <- str_extract(img_pathrows, '[0-9]{3}$')

    if (length(zipfiles) == 0) {
        stop('No images found')
    }

    for (n in 1:length(zipfiles)) {
        zipfile_path <- file.path(in_folder, zipfiles[n])
        # Figure out which satellite the image is from
        year <- format(img_dates[n], '%Y')
        julian_day <- format(img_dates[n], '%j')
        this_out_folder <- file.path(out_folder,
                                     paste0(img_paths[n],'-', img_rows[n], '_', year, 
                                            '-', julian_day, '_', img_sensors[n]))
        if (!file_test('-d', this_out_folder)) {
            dir.create(this_out_folder)
        } else {
            message(paste('Skipping', zipfiles[n], '- output dir', 
                          this_out_folder, 'already exists.'))
            next
        }
        message(paste0(n, ' of ', length(zipfiles), '. Extracting ', zipfiles[n], ' to ', this_out_folder))
        ret_code <- untar(zipfile_path, exdir=file.path(this_out_folder))
        if (ret_code != 0) {
            message(paste('WARNING: error extracting', zipfiles[n], '- return 
                          code', ret_code))
        }
    }
}
