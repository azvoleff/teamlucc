# Function to test if ENVI will load in IDL
check_ENVI_IDL <- function(idl) {
    idl_out <- system(paste(shQuote(idl), '-e "e=ENVI(/HEADLESS)"'), 
                      intern=TRUE)
    if (sum(grepl("Restored file: ENVI", idl_out)) > 0) {
        return(TRUE)
    } else {
        return(FALSE)
    }
}

# Function to ensure only character variables handed to IDL are quoted
format_IDL_param <- function(varname, varvalue) {
    if (is.character(varvalue)) {
        param <- paste0(varname, '="', varvalue, '"\n')
    } else if (is.list(varvalue)) {
        param <- paste0(varname, '=[')
        if (length(varvalue) > 0) {
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
        }
        param <- paste0(param, ']\n')
    } else {
        param <- paste0(varname, '=', varvalue, '\n')
    }
    return(param)
}

#' @importFrom tools file_path_sans_ext
cloud_remove_IDL <- function(cloudy, clear, cloud_mask, out_name,
                             algorithm, num_class, min_pixel, max_pixel, 
                             cloud_nbh, DN_min, DN_max, 
                             verbose, idl, byblock, overwrite,
                             patch_long=1000) {
    if (verbose > 0) {
        warning("verbose not supported with CLOUD_REMOVE and CLOUD_REMOVE_FAST algorithms")
    }
    if (algorithm == 'CLOUD_REMOVE_FAST') {
        script_path <- system.file("idl", "CLOUD_REMOVE_FAST.pro", 
                                   package="teamlucc")
        function_name <- 'CLOUD_REMOVE_FAST'
    } else if (algorithm == 'CLOUD_REMOVE') {
        script_path <- system.file("idl", "CLOUD_REMOVE.pro", 
                                   package="teamlucc")
        function_name <- 'CLOUD_REMOVE'
    } else {
        stop(paste0('unrecognized cloud fill algorithm "', algorithm, '"'))
    }
    
    if (!(file_test('-x', idl) || file_test('-f', idl))) {
        stop('IDL not found - check "idl" parameter')
    }

    if (!check_ENVI_IDL(idl)) {
        stop("Unable to load ENVI in IDL - do you have ENVI and IDL licenses, and ENVI >= 5.0?")
    }

    if (!byblock) {
        patch_long <- max(dim(cloudy)) + 1
    }

    # Save proj4string and extent to ensure the same proj4string and extent is 
    # returned even if they are changed by IDL
    orig_proj <- proj4string(cloudy)
    orig_ext <- extent(cloudy)

    # Write in-memory rasters to files for hand off to IDL. The capture.output 
    # line is used to avoid printing the rasterOptions to screen as they are 
    # temporarily reset.
    dummy <- capture.output(def_format <- rasterOptions()$format)
    rasterOptions(format='ENVI')
    cloudy <- writeRaster(cloudy, rasterTmpFile(), 
                          datatype=dataType(cloudy)[1])
    clear <- writeRaster(clear, rasterTmpFile(), datatype=dataType(clear)[1])
    cloud_mask <- writeRaster(cloud_mask, rasterTmpFile(), 
                              datatype=dataType(cloud_mask)[1])
    cloudy_file <- filename(cloudy)
    clear_file <- filename(clear)
    cloud_mask_file <- filename(cloud_mask)
    dummy <- capture.output(rasterOptions(format=def_format))

    param_names <- c("cloudy_file", "clear_file", "mask_file", "out_name", 
                     "num_class", "min_pixel", "extent1", "DN_min", "DN_max", 
                     "patch_long")
    param_vals <- list(cloudy_file, clear_file, cloud_mask_file, out_name, 
                       num_class, min_pixel, cloud_nbh, DN_min, DN_max, 
                       patch_long)
    idl_params <- mapply(format_IDL_param, param_names, param_vals)
    idl_params <- paste(idl_params, collapse='')

    script_dir <- dirname(script_path)
    idl_script <- tempfile(fileext='.pro')
    idl_cmd <- paste0('CD, "', script_dir, '"\n', idl_params, function_name, ',', 
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

    filled <- brick(out_name)

    filled[filled < DN_min] <- NA

    # Ensure original proj4string and extent are saved and returned
    proj4string(filled) <- orig_proj
    extent(filled) <- orig_ext
    filled <- writeRaster(filled, filename=out_name, overwrite=TRUE, 
                          datatype=dataType(filled)[1])
    return(filled)
}

# Wrapper around C++ cloud fill function, to enable calling the function with 
# rasterEngine
#' @import Rcpp
cloud_fill_rasterengine <- function(cloudy, clear, cloud_mask, algorithm, 
                                    num_class, min_pixel, max_pixel, cloud_nbh, 
                                    DN_min, DN_max, verbose, ...) {
    dims=dim(cloudy)
    # RcppArmadillo crashes when you pass it a cube, so resize and pass 
    # mats
    cloudy <- array(cloudy, dim=c(dims[1] * dims[2], dims[3]))
    clear <- array(clear, dim=c(dims[1] * dims[2], dims[3]))
    cloud_mask <- array(cloud_mask, dim=c(dims[1] * dims[2]))
    filled <- call_cpp_cloud_fill(cloudy, clear, cloud_mask, algorithm, dims, 
                                  num_class,  min_pixel, max_pixel, cloud_nbh, 
                                  DN_min, DN_max, verbose)
    # RcppArmadillo crashes when you return a cube, so resize the returned 
    # mat
    filled <- array(filled, dim=c(dims[1], dims[2], dims[3]))
    return(filled)
}

# This function decides which RcppArmadillo exported function to call: 
# cloud_fill, or cloud_fill_simple
call_cpp_cloud_fill <- function(cloudy, clear, cloud_mask, algorithm, dims, 
                                num_class, min_pixel, max_pixel, cloud_nbh, 
                                DN_min, DN_max, verbose, ...) {
    if (algorithm == "teamlucc") {
        filled <- cloud_fill(cloudy, clear, cloud_mask, dims, num_class, 
                             min_pixel, max_pixel, cloud_nbh, DN_min, DN_max, 
                             verbose)
    } else if (algorithm == "simple") {
        filled <- cloud_fill_simple(cloudy, clear, cloud_mask, dims, num_class, 
                                    cloud_nbh, DN_min, DN_max, verbose)
    } else {
        stop(paste0('unrecognized cloud fill algorithm "', algorithm, '"'))
    }
    return(filled)
}

#' @importFrom spatial.tools rasterEngine
cloud_remove_R <- function(cloudy, clear, cloud_mask, out_name, algorithm, 
                           num_class, min_pixel, max_pixel, cloud_nbh, DN_min, 
                           DN_max, verbose, byblock, overwrite) {
    # Note that call_cpp_cloud_fill uses the algorithm to decide whether to 
    # call cloud_fill or cloud_fill_simple (and call_cpp_cloud_fill is called 
    # by cloud_fill_rasterengine)
    if (byblock) {
        bs <- blockSize(cloudy)
        out <- brick(cloudy, values=FALSE)
        out <- writeStart(out, out_name, overwrite=overwrite)
        for (block_num in 1:bs$n) {
            if (verbose > 0) {
                message("Processing block ", block_num, " of ", bs$n, "...")
            }
            dims <- c(bs$nrows[block_num], ncol(cloudy), nlayers(cloudy))
            cloudy_bl <- array(getValuesBlock(cloudy, row=bs$row[block_num],
                                              nrows=bs$nrows[block_num]),
                               dim=c(dims[1] * dims[2], dims[3]))
            clear_bl <- array(getValuesBlock(clear, row=bs$row[block_num],
                                            nrows=bs$nrows[block_num]),
                            dim=c(dims[1] * dims[2], dims[3]))
            cloud_mask_bl <- array(getValuesBlock(cloud_mask, 
                                                  row=bs$row[block_num], 
                                                  nrows=bs$nrows[block_num]), 
                                   dim=c(dims[1] * dims[2]))
            filled <- call_cpp_cloud_fill(cloudy_bl, clear_bl, cloud_mask_bl, 
                                          algorithm, dims, num_class, 
                                          min_pixel, max_pixel, cloud_nbh, 
                                          DN_min, DN_max, verbose>1)
            out <- writeValues(out, filled, bs$row[block_num])
        }
        out <- writeStop(out)
        # out <- rasterEngine(cloudy=cloudy, clear=clear, 
        # cloud_mask=cloud_mask,
        #                     fun=cloud_fill_rasterengine,
        #                     args=list(algorithm=algorithm, num_class=num_class, 
        #                               min_pixel=min_pixel, max_pixel=max_pixel, 
        #                               cloud_nbh=cloud_nbh, DN_min=DN_min, 
        #                               DN_max=DN_max, verbose=verbose),
        #                     processing_unit='chunk',
        #                     outbands=nlayers(cloudy), outfiles=1,
        #                     verbose=verbose,
        #                     filename=out_name)
    } else {
        dims <- dim(cloudy)
        out_datatype <- dataType(cloudy)[1]
        out <- brick(cloudy, values=FALSE, filename=out_name)
        # RcppArmadillo crashes when you pass it a cube, so resize and pass 
        # mats
        cloudy <- array(getValues(cloudy), dim=c(dims[1] * dims[2], dims[3]))
        clear <- array(getValues(clear), dim=c(dims[1] * dims[2], dims[3]))
        cloud_mask <- array(getValues(cloud_mask), dim=c(dims[1] * dims[2]))
        filled <- call_cpp_cloud_fill(cloudy, clear, cloud_mask, algorithm, 
                                      dims, num_class, min_pixel, max_pixel, 
                                      cloud_nbh, DN_min, DN_max, verbose>1)
        out <- setValues(out, filled)
        out <- writeRaster(out, out_name, datatype=out_datatype, 
                           overwrite=overwrite)
    }

    return(out)
}

#' Remove clouds From Landsat imagery
#'
#' This code uses one of several different algorithms (depending on the 
#' settings, see Details) to fill heavy clouds in a Landsat image.
#'
#' The \code{algorithm} parameter determines what algorithm is used for the 
#' cloud fill. \code{algorithm} must be one of: "CLOUD_REMOVE", 
#' "CLOUD_REMOVE_FAST", "teamlucc", or "simple" (the default). If set to 
#' "CLOUD_REMOVE" the script uses a (slightly modified to be called from R) 
#' version of  Xiaolin Zhu's NSPI IDL code. If set to "CLOUD_REMOVE_FAST", the 
#' algorithm uses the "fast" version of Xiaolin's code. Both of these two 
#' algorithms require an IDL license to run (and therefore \code{idl_path} must 
#' be set).  The "teamlucc" algorithm uses a version of the NSPI algorithm 
#' (based on the CLOUD_REMOVE code) that is coded in C++ and  can be run from R 
#' without an IDL license. The "simple" algorithm uses a cloud fill model that 
#' is based on fitting a linear model to the surface reflectance from the clear 
#' image in a window around each cloud, and using this linear model to predict 
#' reflectance in unobserved (cloudy) areas.
#'
#' @export
#' @param cloudy the cloudy image (base image) as a \code{Raster*}
#' @param clear the clear image as a \code{Raster*} to use for filling 
#' \code{img_cloudy}
#' @param cloud_mask the cloud mask as a \code{RasterLayer}, with each cloud 
#' patch assigned a unique integer code. Areas that are clear in both 
#' \code{cloudy_rast} and \code{clear_rast} should be coded 0, while areas that 
#' are clouded in \code{clear_rast} should be coded -1.
#' @param out_name filename for cloud filled image
#' @param algorithm must be one of: "CLOUD_REMOVE", "CLOUD_REMOVE_FAST", 
#' "teamlucc", or "simple". Default is "simple". See Details.
#' @param num_class set the estimated number of classes in image
#' @param min_pixel the sample size of similar pixels (ignored when 
#' \code{algorithm==TRUE})
#' @param max_pixel the maximum sample size to search for similar pixels 
#' (ignored when \code{algorithm==TRUE})
#' @param cloud_nbh the range of cloud neighborhood (in pixels)
#' @param DN_min the minimum valid DN value (default of 0)
#' @param DN_max the maximum valid DN value (default of 10000 assumes 2 byte 
#' integer imagery)
#' @param idl path to the IDL binary on your machine (on Windows, the path to 
#' idl.exe)
#' @param verbose whether to print detailed status messages. Set to FALSE or 0 
#' for no status messages. Set to 1 for basic status messages. Set to 2 for 
#' detailed status messages.
#' @param byblock whether to process images block by block 
#' (\code{byblock=TRUE}) or all at once (\code{byblock=FALSE}). Use 
#' \code{byblock=FALSE} with caution, as this option will cause the cloud fill 
#' routine to consume a large amount of memory.
#' @param overwrite whether to overwrite \code{out_name} if it already exists
#' @param ... additional arguments passed to the chosen cloud fill routine
#' @return \code{Raster*} with cloud-filled image
#' @references Zhu, X., Gao, F., Liu, D., Chen, J., 2012. A modified
#' neighborhood similar pixel interpolator approach for removing thick clouds 
#' in Landsat images. Geoscience and Remote Sensing Letters, IEEE 9, 521--525.
#' @examples
#' \dontrun{
#' cloudy <- raster(system.file('tests', 'testthat_idl', 'cloud_remove', 
#' 'L20080724_cloudy', package='teamlucc'))
#' clear <- raster(system.file('tests', 'testthat_idl', 'cloud_remove', 
#' 'L20080606', package='teamlucc'))
#' cloud_mask <- raster(system.file('tests', 'testthat_idl', 'cloud_remove', 
#' 'cloud_mask', package='teamlucc'))
#' filled <- cloud_remove(cloudy, clear, cloud_mask, fast=TRUE)
#' }
cloud_remove <- function(cloudy, clear, cloud_mask, out_name=NULL, 
                         algorithm='simple',
                         num_class=4, min_pixel=20, max_pixel=1000, 
                         cloud_nbh=10, DN_min=0, DN_max=10000, 
                         idl="C:/Program Files/Exelis/IDL83/bin/bin.x86_64/idl.exe",
                         verbose=FALSE, byblock=TRUE, overwrite=FALSE, ...) {
    if (!(algorithm %in% c('CLOUD_REMOVE', 'CLOUD_REMOVE_FAST', 'teamlucc', 
                           'simple'))) {
        stop('algorithm must be one of "CLOUD_REMOVE", "CLOUD_REMOVE_FAST", "teamlucc", or "simple"')
    }

    if (verbose > 0) {
        message('Using "', algorithm, '" algorithm.')
    }
    
    if (!(class(cloudy) %in% c("RasterLayer", "RasterStack", "RasterBrick"))) {
        stop('cloudy must be a Raster* object')
    }
    if (!(class(clear) %in% c("RasterLayer", "RasterStack", "RasterBrick"))) {
        stop('clear must be a Raster* object')
    }
    if (!(class(cloud_mask) %in% c("RasterLayer"))) {
        stop('cloud_mask must be a RasterLayer object')
    }
    compareRaster(cloudy, clear)
    if (nlayers(cloudy) != nlayers(clear)) {
        stop('number of layers in cloudy must match number of layers in clear')
    }
    if (nlayers(cloud_mask) != 1) {
        stop('cloud_mask should have only one layer')
    }

    if (is.null(out_name)) {
        out_name <- rasterTmpFile()
    } else {
        out_name <- normalizePath(out_name, mustWork=FALSE)
        if (!file_test('-d', dirname(out_name))) {
            stop('output folder does not exist')
        }
        if (file_test('-f', out_name) & !overwrite) {
            stop('output file already exists - use a different "out_name"')
        }
    }
    
    if (algorithm %in% c('CLOUD_REMOVE', 'CLOUD_REMOVE_FAST')) {
        filled <- cloud_remove_IDL(cloudy, clear, cloud_mask, out_name,
                                   algorithm, num_class, min_pixel, max_pixel, 
                                   cloud_nbh, DN_min, DN_max, verbose, idl, 
                                   byblock, overwrite, ...)
    } else if (algorithm %in% c('teamlucc', 'simple')) {
        filled <- cloud_remove_R(cloudy, clear, cloud_mask, out_name, 
                                 algorithm, num_class, min_pixel, max_pixel, 
                                 cloud_nbh, DN_min, DN_max, verbose, byblock, 
                                 overwrite, ...)
    } else {
        stop(paste0('unrecognized cloud fill algorithm "', algorithm, '"'))
    }

    names(filled) <- names(cloudy)
    return(filled)
}
