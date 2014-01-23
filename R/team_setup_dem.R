#' Setup the DEM mosaic for a particular TEAM site
#'
#' @export
#' @param dem_path a list of digital elevation models (DEMs) that (when 
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
#' @param notify notifier to use (defaults to \code{print} function). See the 
#' \code{notifyR} package for one way of sending notifications from R. The 
#' \code{notify} function should accept a string as the only argument.
#' @examples
#' \dontrun{
#' dem_path <- 'H:/Data/TEAM/VB/Rasters/DEM/ASTER'
#' team_setup_dem(dem_path, "VB", 'H:/Data/TEAM/VB/LCLUC_Analysis/')
#' }
team_setup_dem <- function(dem_path, sitecode, output_path, sample_image=NULL, 
                           n_cpus=1, overwrite=FALSE, notify=print) {
    timer <- Track_time(notify)

    timer <- start_timer(timer, label='Setting up DEMs')

    if (n_cpus > 1) beginCluster(n_cpus)

    # Below should recognize ASTER GDEMV2, SRTM tiles from Earth Explorer, or 
    # CGIAR-SRTM tiles renamed locally with the rename_tiles.R script
    dem_list <- dir(dem_path, pattern='(_3arc_v2\\.bil)|(_dem\\.tif)|(cgiar_srtm_[EW][0-9]{3}_[NS][0-9]{2}.tif$)')
    dem_list <- file.path(dem_path, dem_list)
    dem_rasts <- lapply(dem_list, raster)

    if (length(dem_list) > 1) {
        ######################################################################
        # Verify projections of DEMs match
        dem_prj <- projection(dem_rasts[[1]])
        if (any(lapply(dem_rasts, projection) != dem_prj)) {
            stop("each DEM in dem_list must have the same projection")
        }

        ######################################################################
        # Mosaic DEMs
        timer <- start_timer(timer, label='Mosaicing DEMs')
        # See http://bit.ly/1dJPIeF re issue in raster that necessitates below 
        # workaround
        # TODO: Contact Hijmans re possible fix
        mosaicargs <- dem_rasts
        mosaicargs$fun <- mean
        dem_mosaic <- do.call(mosaic, mosaicargs)
        timer <- stop_timer(timer, label='Mosaicing DEMs')
    } else {
        dem_mosaic <- dem_rasts[[1]]
    }
    dem_mosaic <- round(dem_mosaic)
    dataType(dem_mosaic) <- 'INT2S'

    dem_mosaic_filename <- file.path(output_path,
                                     paste0(sitecode, '_dem_mosaic.envi'))
    if (is.null(sample_image) | compareRaster(sample_image, dem_mosaic,
                                              extent=FALSE, rowcol=FALSE, 
                                              crs=TRUE, res=TRUE, orig=TRUE, 
                                              stopiffalse=FALSE)) {
        dem_mosaic <- writeRaster(dem_mosaic, dem_mosaic_filename, 
                                  overwrite=overwrite, 
                                  datatype=dataType(dem_mosaic))
    } else {
        timer <- start_timer(timer, label='Reprojecting DEM mosaic')
        # The below lines construct to_ext as the extent the image will be 
        # projected to. This extent must cover the same area as the dem_mosaic, 
        # but must have the same resolution, CRS and origin as the sample_image
        to_ext <- projectExtent(dem_mosaic, crs(sample_image))
        to_res <- res(sample_image)
        xmin(to_ext) <- floor(xmin(to_ext) / to_res[1]) * to_res[1]
        xmax(to_ext) <- ceiling(xmax(to_ext) / to_res[1]) * to_res[1]
        ymin(to_ext) <- floor(ymin(to_ext) / to_res[2]) * to_res[2]
        ymax(to_ext) <- ceiling(ymax(to_ext) / to_res[2]) * to_res[2]
        res(to_ext) <- to_res
        dem_mosaic <- projectRaster(from=dem_mosaic, to=to_ext,
                                    filename=dem_mosaic_filename, 
                                    overwrite=overwrite, 
                                    datatype=dataType(dem_mosaic))
        timer <- stop_timer(timer, label='Reprojecting DEM mosaic')
    }

    timer <- start_timer(timer, label='Calculating slope and aspect')
    slopeaspect_filename <- file.path(output_path,
                                     paste0(sitecode, '_dem_mosaic_slopeaspect.envi'))
    # Note that the default output of 'terrain' is in radians
    slopeaspect <- terrain(dem_mosaic, opt=c('slope', 'aspect'))
    slopeaspect$aspect <- calc(slopeaspect$aspect, fun=function(vals) {
        vals[vals >= 2*pi] <- 0
        vals
        })
    # Note that slopeaspect is scaled - slope by 10000, and aspect by 1000 so 
    # that the layers can be saved as INT2S
    slopeaspect <- stack(round(raster(slopeaspect, layer=1) * 10000),
                         round(raster(slopeaspect, layer=2) * 1000))
    slopeaspect <- writeRaster(slopeaspect, filename=slopeaspect_filename, 
                               overwrite=overwrite, datatype='INT2S')
    timer <- stop_timer(timer, label='Calculating slope and aspect')

    if (n_cpus > 1) endCluster()

    timer <- stop_timer(timer, label='Setting up DEMs')
}
