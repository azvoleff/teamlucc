#' Normalize a set of preprocessed CDR images to a base image
#'
#' This function uses model II regression to perform relative normalization to 
#' match a set of Landsat CDR surface reflectance images. The function assumes 
#' the images were preprocessed using the \code{auto_preprocess_landsat} 
#' function. A base image can be optionally supplied. If a base image is not 
#' supplied, then the function will calculate the percent cloud cover of each 
#' input image, and automatically choose the image with the least cloud cover 
#' as the base image. This function assumes that the images (and image masks) 
#' were preprocessed using the \code{auto_preprocess_landsat} function.
#'
#' This function will run in parallel if a parallel backend is registered with 
#' \code{\link{foreach}}.
#'
#' @export
#' @import foreach
#' @importFrom tools file_path_sans_ext
#' @param image_files list of filenames for images to normalize
#' @param base (optional) filename of base image. If not supplied, the base 
#' image will be automatically chosen from the images in \code{image_files}, as 
#' the image with the lowest percent cloud cover.
#' @param overwrite whether to overwrite existing files
#' @return nothing - used for side effect of normalizing imagery
auto_normalize <- function(image_files, base, overwrite=FALSE) {
    stopifnot(length(image_files) >= 1)

    image_stacks <- lapply(image_files, stack)

    mask_files <- paste0(file_path_sans_ext(image_files), '_masks', 
                         extension(image_files))
    mask_stacks <- lapply(mask_files, stack)

    if (!missing(base)) {
        base_img_file <- base
    } else if (missing(base) & (length(image_files) == 1)) {
        stop('length of image_files is 1 but no base image was supplied')
    } else {
        # Figure out which image has lowest percent cloud cover - use that 
        # image as the base image
        pct_clouds <- function(cloud_mask) {
            clouded_pixels <- calc(cloud_mask, fun=function(vals) {
                # For fmask layer, 2 is cloud shadow, and 4 is cloud
                (vals == 2) | (vals == 4)
            })
            num_clouds <- cellStats(clouded_pixels, stat='sum', na.rm=TRUE)
            # For fmask layer, 255 is fill
            num_clear <- cellStats(cloud_mask != 255, stat='sum', na.rm=TRUE)
            return((num_clouds / (num_clouds + num_clear)) * 100)
        }

        cloud_cover <- foreach(mask_stack=iter(mask_stacks),
                 .packages=c('teamlucc', 'stringr', 'rgdal'),
                 .combine=c) %dopar% {
            # Note that fmask layer is 2nd layer in stack
            pct_clouds(mask_stack[[2]])
        }
        base_index <- which(cloud_cover == min(cloud_cover))
        
        base_img <- image_stacks[[base_index]]
        image_stacks <- image_stacks[-base_index]

        base_img_file <- image_files[[base_index]]
        image_files <- image_files[-base_index]

        base_mask <- mask_stacks[[base_index]]
        mask_stacks <- mask_stacks[-base_index]
    }

    # Copy the base image to a new file with the _base.tif extension
    base_copy_filename <- paste0(file_path_sans_ext(base_img_file), 
                                 '_normbase.tif')
    base_img <- writeRaster(base_img, filename=base_copy_filename, 
                            datatype=dataType(base_img)[1], 
                            overwrite=overwrite)
    base_mask_copy_filename <- paste0(file_path_sans_ext(base_img_file), 
                                      '_normbase_masks.tif')
    base_mask <- writeRaster(base_mask, filename=base_mask_copy_filename, 
                             datatype=dataType(base_mask)[1], 
                             overwrite=overwrite)

    stopifnot(length(image_files) == length(image_stacks))
    stopifnot(length(image_files) == length(mask_stacks))

    image_file=image_stack=NULL
    # Now normalize each remaining image to the _base.tif file
    foreach (image_file=iter(image_files), image_stack=iter(image_stacks), 
             mask_stack=iter(mask_stacks),
             .packages=c('teamlucc', 'stringr', 'tools')) %dopar% {
        message(paste('Preprocessing ', image_file))
        output_normed_file <- paste0(file_path_sans_ext(image_file), 
                                     '_normalized.tif')
        output_normed_masks_file <- paste0(file_path_sans_ext(image_file), 
                                           '_normalized_masks.tif')

        # Note that fmask layer is 2nd layer in stack
        missing_vals <- overlay(base_mask[[2]], mask_stack[[2]],
                            fun=function(base_vals, this_vals) {
            # Only use clear pixels when normalizing (0 in fmask)
            (base_vals != 0) & (this_vals != 0)
        }, datatype=dataType(base_mask))

        if (ncell(image_stack) > 500000) {
            size <- 500000
        } else {
            size <- ncell(image_stack)
        }
        normed_image <- normalize(base_img, image_stack, missing_vals, size=size)

        normed_image <- writeRaster(normed_image, filename=output_normed_file, 
                                    datatype=dataType(base_img)[1], 
                                    overwrite=overwrite)
        mask_stack <- writeRaster(mask_stack, 
                                  filename=output_normed_masks_file, 
                                  datatype=dataType(mask_stack)[1], 
                                  overwrite=overwrite)
    }
}
