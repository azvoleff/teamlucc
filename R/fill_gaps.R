#' Perform SLC-off gap fill of Landsat 7 ETM+ image
#'
#' Calls GNSPI.pro IDL script by Xiaolin Zhu to fill gaps in SLC-off Landsat 7 
#' ETM+ image. The script requires \code{fill} to be a TM (Landsat 5) image.  
#' \code{slc_off} must be a Landsat 7 SLC-off image.
#'
#' If supplied, \code{timeseries} should be a list of TM images.  Performing 
#' gap fill using SLC-Off ETM+ images as the input is not yet supported.
#' 
#' Pixels in gaps, background, and/or clouds in \code{slc_off}, 
#' \code{input_image}, and the images in \code{timeseries} should be coded as 
#' 0.
#'
#' @importFrom tools file_path_sans_ext
#' @param slc_off the SLC-off Landsat 7 file to gap fill, as a \code{Raster*}
#' @param fill the first TM image to use to fill in the gaps, as a 
#' \code{Raster*}
#' @param timeseries a timeseries of TM images as \code{Raster*} objects to use 
#' as additional inputs to the gap fill algorithm (optional)
#' @param out_base path and base filename for the output file. The script will 
#' save the output files by appending "_GNSPI.envi" and 
#' "_GNSPI_uncertainty.envi" to this base filename.
#' @param ext file extension to use when and when saving output rasters 
#' (determines output file format). Must be supported by 
#' \code{\link{writeRaster}}.
#' @param algorithm the algorithm to use, as a string ("GNSPI_IDL" is currently 
#' the only supported algorithm)
#' @param sample_size the sample size of sample pixels
#' @param size_wind the maximum window size
#' @param class_num the estimated number of classes
#' @param DN_min the minimum DN value of the image
#' @param DN_max the maximum DN value of the image
#' @param patch_long the size of block, to process whole ETM scene, set to 1000
#' @param idl path to the IDL binary
#' @param verbose whether to print detailed status messages
#' @param overwrite whether to overwrite output files if they already exist
#' @return a list of two rasters: 1) "filled", the gap filled image, and 2) 
#' "uncertainty", the uncertainty image.
#' @export
#' @references Zhu, X., Liu, D., Chen, J., 2012. A new geostatistical approach 
#' for filling gaps in Landsat ETM+ SLC-off images. Remote Sensing of 
#' Environment 124, 49--60.
#' @examples
#' \dontrun{
#' slc_off <- brick(system.file('tests', 'testthat_idl', 'fill_gaps', 
#' 'TM20100429_toaR_gap', package='teamlucc'))
#' fill <- brick(system.file('tests', 'testthat_idl', 'fill_gaps', 
#' 'TM20100515_toaR', package='teamlucc'))
#' timeseries <- c(brick(system.file('tests', 'testthat_idl', 'fill_gaps', 
#' 'TM20100208_toaR', package='teamlucc')))
#' filled <- fill_gaps(slc_off, fill, timeseries)
#' }
fill_gaps <- function(slc_off, fill, timeseries=c(), out_base=NULL, ext="tif",
                      algorithm="GNSPI_IDL", sample_size=20, size_wind=12, 
                      class_num=4, DN_min=0.0, 
                      DN_max=1.0, patch_long=1000,
                      idl="C:/Program Files/Exelis/IDL83/bin/bin.x86_64/idl.exe",
                      verbose=FALSE, overwrite=FALSE) {
    if (!(class(slc_off) %in% c("RasterLayer", "RasterStack", "RasterBrick"))) {
        stop('slc_off must be a Raster* object')
    }
    if (!(class(fill) %in% c("RasterLayer", "RasterStack", "RasterBrick"))) {
        stop('fill must be a Raster* object')
    }
    if (nlayers(slc_off) != nlayers(fill)) {
        stop('number of layers in slc_off must match number of layers in fill')
    }
    compareRaster(slc_off, fill)
    if (!(algorithm %in% c('GNSPI_IDL'))) {
        stop('algorithm must be "GNSPI_IDL" - no other algorithms are supported')
    }

    ext <- gsub('^[.]', '', ext)

    for (timeseries_img in timeseries) {
        if (!(class(timeseries_img) %in% c("RasterLayer", "RasterStack", "RasterBrick"))) {
            stop('each timeseries image be a Raster* object')
        }
        if (nlayers(slc_off) != nlayers(timeseries_img)) {
            stop('number of layers in slc_off must match number of layers of each image in timeseries')
        }
        compareRaster(slc_off, timeseries_img)
    }

    if (is.null(out_base)) {
        out_base <- file_path_sans_ext(rasterTmpFile())
    } else {
        out_base <- normalizePath(out_base, mustWork=FALSE)
        if (!file_test('-d', dirname(out_base))) {
            stop('output folder does not exist')
        }
        if (!overwrite && file_test('-f', file.path(out_base, paste0('_GNSPI.', ext)))) {
            stop('output file already exists - use a different "out_base"')
        }
        if (!overwrite && file_test('-f', file.path(out_base, paste0('_GNSPI_uncertainty.', ext)))) {
            stop('output uncertainty file already exists - use a different "out_base"')
        }
    }
    
    if (algorithm == "GNSPI_IDL") {
        filled <- fill_gaps_idl(slc_off, fill, timeseries, out_base, 
                                sample_size, size_wind, class_num, DN_min, 
                                DN_max, patch_long, idl, algorithm, ext, 
                                verbose)
    } else {
        stop("Native R gap filling not yet supported")
    }

    return(filled)
}

fill_gaps_idl <- function(slc_off, fill, timeseries, out_base, sample_size, 
                          size_wind, class_num, DN_min, DN_max, patch_long, 
                          idl, algorithm, ext, verbose) {
    if (verbose) {
        warning('verbose=TRUE not supported when algorithm="GNSPI_IDL"')
    }

    script_path <- system.file("idl", "GNSPI.pro", package="teamlucc")
    if (!(file_test('-x', idl) || file_test('-f', idl))) {
        stop('IDL not found - check "idl" parameter')
    }

    if (!check_ENVI_IDL(idl)) {
        stop("Unable to load ENVI in IDL - do you have ENVI and IDL licenses, and ENVI >= 5.0?")
    }

    # Save proj4string and extend to ensure the same proj4string and extent is 
    # returned even if they are changed by IDL
    orig_proj <- proj4string(slc_off)
    orig_ext <- extent(slc_off)
    orig_datatype <- dataType(slc_off)[1]

    # Write in-memory rasters to files for hand off to IDL. The capture.output 
    # line is used to avoid printing the rasterOptions to screen as they are 
    # temporarily reset.
    dummy <- capture.output(def_format <- rasterOptions()$format)
    rasterOptions(format='ENVI')
    slc_off <- writeRaster(slc_off, rasterTmpFile(), datatype=dataType(slc_off)[1])
    slc_off_file <- filename(slc_off)
    fill <- writeRaster(fill, rasterTmpFile(), datatype=dataType(fill)[1])
    fill_file <- filename(fill)
    timeseries_files <- c()
    if (length(timeseries) == 0) {
        # Ensure timeseries_files equals an empty matrix, in IDL format
        timeseries_files <- list()
    } else {
        for (timeseries_img in timeseries) {
            timeseries_img <- writeRaster(timeseries_img, rasterTmpFile(), datatype=dataType(fill)[1])
            timeseries_files <- c(timeseries_files, filename(timeseries_img))
        }
    }
    temp_dir <- tempdir()
    dummy <- capture.output(rasterOptions(format=def_format))

    # Save IDL output to a temp folder - it will be copied over and saved with 
    # writeRaster later to ensure the extents and projection are not modified 
    # from those of the original files.
    temp_out_base <- file_path_sans_ext(rasterTmpFile())
    param_vals <- list(slc_off_file, fill_file, timeseries_files,
                       temp_out_base, sample_size, size_wind, class_num,
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

    filled <- brick(paste0(temp_out_base, '_GNSPI.envi'))
    filled_out_file <- paste0(out_base, paste0('_GNSPI.', ext))
    proj4string(filled) <- orig_proj
    extent(filled) <- orig_ext
    filled <- writeRaster(filled, filename=filled_out_file, overwrite=TRUE, 
                          datatype=orig_datatype)

    uncertainty <- brick(paste0(temp_out_base, '_GNSPI_uncertainty.envi'))
    uncertainty_out_file <- paste0(out_base, paste0('_GNSPI_uncertainty.', ext))
    proj4string(uncertainty) <- orig_proj
    extent(uncertainty) <- orig_ext
    uncertainty <- writeRaster(uncertainty, filename=uncertainty_out_file, 
                               overwrite=TRUE, datatype=orig_datatype)

    return(list(filled=filled, uncertainty=uncertainty))
}
