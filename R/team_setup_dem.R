#' Setup the DEM mosaic for a particular TEAM site
#'
#' @export
#' @import spatial.tools
#' @param dem_list a list of digital elevation models (DEMs) that (when 
#' mosaiced) covers the full extent of all the images in the image_list.
#' @param sitecode code to use as a prefix for all filenames
#' @param output_path the path to use for the output 
#' @param sample_image a \code{Raster*} or the path to an image in a format 
#' readable by \code{raster}. The sample image will be used to set the 
#' projection system for the output DEM mosaic.
#' @param n_cpus the number of CPUs to use for processes that can run in 
#' parallel
#' @param overwrite whether to overwrite existing files (otherwise an error 
#' will be raised)
#' @examples
#' \dontrun{
#' dem_list <- list('H:/Data/TEAM/VB/Rasters/DEM/ASTER/ASTGTM2_N09W084_dem.tif',
#'                  'H:/Data/TEAM/VB/Rasters/DEM/ASTER/ASTGTM2_N09W085_dem.tif',
#'                  'H:/Data/TEAM/VB/Rasters/DEM/ASTER/ASTGTM2_N10W084_dem.tif',
#'                  'H:/Data/TEAM/VB/Rasters/DEM/ASTER/ASTGTM2_N09W085_dem.tif')
#' team_setup_dem(dem_list, "VB", 'H:/Data/TEAM/VB/LCLUC_Analysis/')
#' }
team_setup_dem <- function(dem_list, sitecode, output_path, sample_image=NULL, 
                            n_cpus=1, overwrite=FALSE) {
    if (n_cpus > 1) sfQuickInit(n_cpus)
    ################################################################################
    # Verify projections of DEMs match
    dem_rasts <- lapply(dem_list, raster)

    dem_prj <- projection(dem_rasts[[1]])
    if (any(lapply(dem_rasts, projection) != dem_prj)) {
        stop("each DEM in dem_list must have the same projection")
    }

    ################################################################################
    # Mosaic DEMs
    print('Mosaicing DEMs...')
    trackTime(action='start')
    # See http://bit.ly/1dJPIeF re issue in raster that necessitates below workaround
    # TODO: Contact Hijmans re possible fix
    mosaicargs <- dem_rasts
    mosaicargs$fun <- mean
    dem_mosaic <- do.call(mosaic, mosaicargs)
    dem_mosaic_filename <- file.path(output_path, paste0(sitecode, '_dem_mosaic.envi'))
    sample_image <- raster(sample_image)
    if (is.null(sample_image) | (projection(sample_image) == projection(dem_mosaic))) {
        dem_mosaic <- writeRaster(dem_mosaic, dem_mosaic_filename, overwrite=overwriteoverwrite)
    } else {
        dem_mosaic <- projectRaster(dem_mosaic, sample_image, filename=dem_mosaic_filename, overwrite=overwrite)
    }
    trackTime()

    print('Running slopeasp_seq...')
    trackTime(action='start')
    slopeaspect_filename <- file.path(output_path, paste0(sitecode, '_dem_mosaic_slopeaspect.envi'))
    slopeaspect <- slopeasp_seq(dem_mosaic, filename=slopeaspect_filename, overwrite=overwrite)
    trackTime()
    if (n_cpus > 1) sfQuickStop()
}
