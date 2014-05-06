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

# Function to prepare fmask band as a cloud mask
prep_fmask <- function(image_dir) {
    lndsr_regex <- '^lndsr.((LT4)|(LT5)|(LE7)|(LE8))[0-9]{6}[12][0-9]{6}[a-zA-Z]{3}[0-9]{2}'
    lndsr_files <- file.path(image_dir, dir(image_dir, pattern=lndsr_regex))
    fmask_band <- raster(lndsr_files[grepl('fmask_band.envi$', lndsr_files)])

    clouds <- (fmask_band == 2) | (fmask_band == 4)
}

#' @importFrom tools file_path_sans_ext
cloud_remove_IDL <- function(cloudy, clear, cloud_mask, out_name,
                             fast, num_class, min_pixel, max_pixel, 
                             cloud_nbh, DN_min, DN_max, 
                             verbose, idl, patch_long=1000, ...) {
    if (verbose) {
        warning("verbose=TRUE not supported when use_IDL=TRUE")
    }
    if (fast) {
        script_path <- system.file("idl", "CLOUD_REMOVE_FAST.pro", 
                                   package="teamlucc")
        function_name <- 'CLOUD_REMOVE_FAST'
    } else {
        script_path <- system.file("idl", "CLOUD_REMOVE.pro", 
                                   package="teamlucc")
        function_name <- 'CLOUD_REMOVE'
    }
    
    if (!(file_test('-x', idl) || file_test('-f', idl))) {
        stop('IDL not found - check "idl" parameter')
    }

    # Write in-memory rasters to files for hand off to IDL. The capture.output 
    # line is used to avoid printing the rasterOptions to screen as they are 
    # temporarily reset.
    dummy <- capture.output(def_format <- rasterOptions()$format)
    rasterOptions(format='ENVI')
    cloudy <- writeRaster(cloudy, rasterTmpFile(), datatype=dataType(cloudy)[1])
    clear <- writeRaster(clear, rasterTmpFile(), datatype=dataType(clear)[1])
    cloud_mask <- writeRaster(cloud_mask, rasterTmpFile(), datatype=dataType(cloud_mask)[1])
    cloudy_file <- filename(cloudy)
    clear_file <- filename(clear)
    cloud_mask_file <- filename(cloud_mask)
    if (is.null(out_name)) {
        out_name <- rasterTmpFile()
    } else {
        out_name <- normalizePath(out_name, mustWork=FALSE)
    }
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

    return(brick(out_name))
}

# Wrapper around C++ cloud fill function, to enable calling the function with 
# rasterEngine
#' @import Rcpp
cloud_fill_cpp_wrapper <- function(cloudy, clear, cloud_mask, num_class, 
                                   min_pixel, max_pixel, cloud_nbh, DN_min, 
                                   DN_max, verbose, ...) {
    dims=dim(cloudy)
    # RcppArmadillo crashes when you pass it a cube, so resize and pass 
    # mats
    cloudy <- array(cloudy, dim=c(dims[1] * dims[2], dims[3]))
    clear <- array(clear, dim=c(dims[1] * dims[2], dims[3]))
    cloud_mask <- array(cloud_mask, dim=c(dims[1] * dims[2]))
    filled <- cloud_fill(cloudy, clear, cloud_mask, dims, num_class, 
                         min_pixel, max_pixel, cloud_nbh, DN_min, DN_max, 
                         verbose)
    # RcppArmadillo crashes when you return a cube, so resize the returned 
    # mat
    filled <- array(filled, dim=c(dims[1], dims[2], dims[3]))
}

#' @importFrom spatial.tools rasterEngine
cloud_remove_R <- function(cloudy, clear, cloud_mask, out_name, fast, 
                           num_class, min_pixel, max_pixel, cloud_nbh, DN_min, 
                           DN_max, verbose, ...) {
    if (fast) {
        stop("fast=TRUE not yet supported when use_IDL=FALSE")
    }

    warning("*** use_IDL=FALSE is still experimental - use results with caution ***")

    # bs <- blockSize(cloudy)
    # out <- brick(cloudy, values=FALSE)
    # out <- writeStart(out, rasterTmpFile())
    # message(paste0("**", bs$n, " blocks to process**"))
    # for (block_num in 1:bs$n) {
    #     message("Processing block ", block_num, " - ", appendLF=FALSE)
    #     cloudy_bl <- getValuesBlock(cloudy, row=bs$row[block_num], nrows=bs$nrows[block_num])
    #     clear_bl <- getValuesBlock(clear, row=bs$row[block_num], nrows=bs$nrows[block_num])
    #     cloud_mask_bl <- getValuesBlock(cloud_mask, row=bs$row[block_num], nrows=bs$nrows[block_num])

    #     dims <- c(bs$nrows[block_num], ncol(cloudy), nlayers(cloudy))

    #     filled <- cloud_fill(cloudy_bl, clear_bl, cloud_mask_bl, dims, 
    #                          num_class, min_pixel, max_pixel, cloud_nbh, 
    #                          DN_min, DN_max)

    #     out <- writeValues(out, filled, bs$row[block_num])
    # }
    # out <- writeStop(out)
    
    if (!is.null(out_name)) {
        out_name <- normalizePath(out_name, mustWork=FALSE)
    }

    out <- rasterEngine(cloudy=cloudy, clear=clear, 
                        cloud_mask=cloud_mask,
                        fun=cloud_fill_cpp_wrapper,
                        args=list(num_class=num_class, min_pixel=min_pixel, 
                        max_pixel=max_pixel, cloud_nbh=cloud_nbh, 
                        DN_min=DN_min, DN_max=DN_max, verbose=verbose),
                        processing_unit='chunk',
                        outbands=nlayers(cloudy), outfiles=1,
                        setMinMax=TRUE,
                        filename=out_name, ...)

    return(out)
}

#' Remove clouds using Xiaolin Zhu's NSPI algorithm
#'
#' This code uses the Neighborhood Similar Pixel Interpolator (NSPI) algorithm 
#' by Xiaolin Zhu to fill heavy clouds in a Landsat image.
#'
#' This code can use either a Xiaolin Zhu's original IDL code, or an R/C++ 
#' implementation native to \code{teamlucc}. The \code{use_IDL} parameter 
#' (defaults to TRUE) decides whether to use Xiaolin's code 
#' (\code{use_IDL=TRUE}) or the  or the R and C++ implementation native to 
#' \code{teamlucc} (\code{use_IDL=FALSE}).
#'
#' The results from running the two alternative versions of the cloud removal 
#' algorithm should be identical. However, there is one important difference to 
#' note: the R/C++ implementation of cloud removal does not currently implement 
#' the \code{fast=TRUE} option.
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
#' @param use_IDL whether to use Xiaolin Zhu's original IDL scripts 
#' (\code{use_IDL=TRUE}) or the C++ implementation native to \code{teamlucc} 
#' (\code{use_IDL=FALSE})
#' @param fast if \code{TRUE}, use the CLOUD_REMOVE_FAST.pro script. If 
#' \code{FALSE}, use the CLOUD_REMOVE.pro script.
#' @param num_class set the estimated number of classes in image
#' @param min_pixel the sample size of similar pixels
#' @param max_pixel the maximum sample size to search for similar pixels
#' @param cloud_nbh the range of cloud neighborhood (in pixels)
#' @param DN_min the minimum valid DN value
#' @param DN_max the maximum valid DN value
#' @param idl path to the IDL binary on your machine (on Windows, the path to 
#' idl.exe)
#' @param verbose whether to print detailed status messages
#' @param ... additional arguments to pass to \code{rasterEngine}
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
                         use_IDL=TRUE, fast=FALSE, num_class=4, min_pixel=20, 
                         max_pixel=1000, cloud_nbh=1, DN_min=0, DN_max=255, 
                         idl="C:/Program Files/Exelis/IDL83/bin/bin.x86_64/idl.exe",
                         verbose=FALSE) {
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
    
    if (use_IDL) {
        filled <- cloud_remove_IDL(cloudy, clear, cloud_mask, out_name,
                                   fast, num_class, min_pixel, max_pixel, 
                                   cloud_nbh, DN_min, DN_max, verbose, idl)
    } else {
        filled <- cloud_remove_R(cloudy, clear, cloud_mask, out_name,
                                 fast, num_class, min_pixel, max_pixel, 
                                 cloud_nbh, DN_min, DN_max, verbose)
    }

    return(filled)
}
