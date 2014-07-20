normalize_extent <- function(te, res=c(30, 30)) {
    # Setup xmin
    te[1] <- round(te[1] - te[1] %% res[1])
    # Setup ymin
    te[2] <- round(te[2] - te[2] %% res[2])
    # Setup xmax
    te[3] <- round(te[3] + res[1] - te[3] %% res[1])
    # Setup ymax
    te[4] <- round(te[4] + res[2] - te[4] %% res[2])
    stopifnot(round(te[1] / res[1]) == (te[1] / res[1]))
    stopifnot(round(te[2] / res[2]) == (te[2] / res[2]))
    stopifnot(round(te[3] / res[1]) == (te[3] / res[1]))
    stopifnot(round(te[4] / res[2]) == (te[4] / res[2]))
    return(te)
}

#' Setup the DEM mosaic for a given AOI
#'
#' This function will setup a set of DEM tiles for each the Landsat path/row 
#' needed to cover a given AOI. The tiles can optionally be cropped to cover 
#' only the portion of each path/row that is included in the AOI, or can cover 
#' the full scene for each path/row needed to cover the AOI.
#'
#' This function uses \code{gdalUtils}, which requires a local GDAL 
#' installation.  See http://trac.osgeo.org/gdal/wiki/DownloadingGdalBinaries 
#' or http://trac.osgeo.org/osgeo4w/ to download the appropriate installer for 
#' your operating system.
#'
#' @export
#' @importFrom wrspathrow pathrow_num
#' @importFrom rgdal readOGR writeOGR
#' @importFrom sp spTransform is.projected
#' @importFrom rgeos gBuffer gIntersects gUnaryUnion gIntersection
#' @importFrom tools file_path_sans_ext
#' @importFrom gdalUtils mosaic_rasters gdalwarp
#' @param aoi area of interest (AOI), as a \code{SpatialPolygonsDataFrame}, to 
#' use as as bounding box when selecting DEMs. Also used to crop and set 
#' projection of the output DEM(s) if \code{crop_to_aoi=TRUE}. Must be in a 
#' projected coordinate system.
#' @param output_path the path to use for the output
#' @param dem_extents a \code{SpatialPolygonsDataFrame} of the extents and 
#' filenames for a set of locally available DEM raster(s) that cover the 
#' \code{aoi}. See the \code{\link{get_extent_polys}} function for one means of 
#' generating this list. \code{dem_extents} must have a "filename" column.
#' @param of output format to use when saving output rasters. See description 
#' of \code{of} in \code{\link{gdalwarp}}.
#' @param ext file extension to use when saving output rasters (determines 
#' output file format). Should match file extension for output format chosen by 
#' \code{of}.
#' @param n_cpus the number of CPUs to use for processes that can run in 
#' parallel
#' @param overwrite whether to overwrite existing files (otherwise an error 
#' will be raised)
#' @param crop_to_aoi whether to crop the dem to the supplied AOI, or to the
#' Landsat path/row polygon for that particular path/row
#' @param notify notifier to use (defaults to \code{print} function). See the 
#' \code{notifyR} package for one way of sending notifications from R. The 
#' \code{notify} function should accept a string as the only argument.
#' @param verbose whether to print detailed status messages and timing 
#' information
#' @return nothing - used for the side effect of setting up DEMs
auto_setup_dem <- function(aoi, output_path, dem_extents, of="GTiff", 
                           ext='tif', n_cpus=1, overwrite=FALSE, 
                           crop_to_aoi=FALSE, notify=print, verbose=FALSE) {
    if (!file_test("-d", output_path)) {
        stop(paste(output_path, "does not exist"))
    }

    if (length(aoi) > 1) {
        stop('aoi should be a SpatialPolygonsDataFrame of length 1')
    }
    stopifnot(is.projected(aoi))

    ext <- gsub('^[.]', '', ext)

    timer <- Track_time(notify)

    pathrows <- pathrow_num(aoi, wrs_type=2, wrs_mode='D', as_polys=TRUE)
    aoi_prproj <- spTransform(aoi, CRS(proj4string(pathrows)))

    timer <- start_timer(timer, label=paste('Processing DEMS for', nrow(pathrows), 
                                            'path/rows'))
    if (crop_to_aoi) {
        # Do a rough crop of the pathrows to the AOI in the pathrow CRS 
        # (pathrow will later be cropped in the AOI CRS).
        pathrows_cropped <- gIntersection(pathrows, aoi_prproj, byid=TRUE)
        row.names(pathrows_cropped) <- row.names(pathrows)
        pathrows_cropped <- SpatialPolygonsDataFrame(pathrows_cropped, 
                                                     data=pathrows@data)
    } else {
        pathrows_cropped <- pathrows
    }

    # Add a 500 m buffer in UTM coordinate system, as 1) slope calculation 
    # requires a window of pixels, and 2) this buffer also helps avoid missing 
    # pixels on the sides of the DEM due to slight misalignments from the 
    # reprojection that will occur later. After buffering transform back to 
    # WGS84 to use for preliminary cropping of the dem mosaic. 
    pathrows_utm <- spTransform(pathrows_cropped,
                                CRS(utm_zone(pathrows_cropped, proj4string=TRUE)))
    pathrows_buffered <- spTransform(gBuffer(pathrows_utm, width=500, byid=TRUE), 
                                 CRS(proj4string(dem_extents)))
    intersecting <- as.logical(gIntersects(dem_extents, 
                                           gUnaryUnion(pathrows_buffered), byid=TRUE))
    if (sum(intersecting) == 0) {
        stop('no intersecting dem extents found')
    } else {
        dem_extents <- dem_extents[intersecting, ]
    }

    dem_list <- dem_extents$filename
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
        if (verbose) timer <- start_timer(timer, label='Mosaicking DEMs')
        mosaic_file <- extension(rasterTmpFile(), ext)
        # Calculate minimum bounding box coordinates:
        mosaic_te <- as.numeric(bbox(pathrows_buffered))
        # Use mosaic_rasters from gdalUtils for speed:
        mosaic_rasters(dem_list, mosaic_file, te=mosaic_te, of=of, 
                       overwrite=overwrite, ot='Int16')
        dem_mosaic <- raster(mosaic_file)
        if (verbose) timer <- stop_timer(timer, label='Mosaicking DEMs')
    } else {
        dem_mosaic <- dem_rasts[[1]]
        mosaic_file <- filename(dem_mosaic)
    }

    for (n in 1:length(pathrows)) {
        pathrow <- pathrows[n, ]
        pathrow_label <- paste(sprintf('%03i', pathrow@data$PATH), 
                               sprintf('%03i', pathrow@data$ROW), sep='-')
        timer <- start_timer(timer, label=paste0('Processing ', n, ' of ', 
                                                 nrow(pathrows), ': ', 
                                                 pathrow_label))

        if (verbose) timer <- start_timer(timer,
                                          label=paste('Cropping/reprojecting DEM mosaic crop for', 
                                          pathrow_label))
        if (crop_to_aoi) {
            to_srs <- proj4string(aoi)
            pathrow_tosrs <- spTransform(pathrow, CRS(to_srs))
            to_ext <- extent(gIntersection(pathrow_tosrs, aoi, byid=TRUE))
        } else {
            to_srs <- utm_zone(pathrow, proj4string=TRUE)
            to_ext <- projectExtent(pathrow, to_srs)
        }
        dem_te <- as.numeric(bbox(to_ext))

        # Ensure origin is set at 0,0
        to_res <- c(30, 30)
        dem_te <- normalize_extent(dem_te, to_res)

        # Calculate minimum bounding box coordinates:
        dem_mosaic_crop_filename <- file.path(output_path,
                                         paste0('dem_', pathrow_label, 
                                                '.', ext))
        dem_mosaic_crop <- gdalwarp(mosaic_file, 
                                    dstfile=dem_mosaic_crop_filename,
                                    te=dem_te, t_srs=to_srs, tr=to_res, 
                                    r='cubicspline', output_Raster=TRUE, 
                                    multi=TRUE, of=of,
                                    wo=paste0("NUM_THREADS=", n_cpus), 
                                    overwrite=overwrite, ot='Int16')
        if (verbose) timer <- stop_timer(timer,
                                         label=paste('Cropping/reprojecting DEM mosaic crop for', 
                                         pathrow_label))

        if (verbose) timer <- start_timer(timer, label=paste('Calculating slope/aspect for', 
                                                pathrow_label))
        slopeaspect_filename <- file.path(output_path,
                                          paste0('slopeaspect_',
                                                 pathrow_label, '.', ext))
        # Note that the default output of 'terrain' is in radians
        slopeaspect <- terrain(dem_mosaic_crop, opt=c('slope', 'aspect'))
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
        if (verbose) timer <- stop_timer(timer, label=paste('Calculating slope/aspect for', 
                                                pathrow_label))
        timer <- stop_timer(timer, label=paste0('Processing ', n, ' of ', 
                                                nrow(pathrows), ': ', 
                                                pathrow_label))
    }

    timer <- stop_timer(timer, label=paste('Processing DEMS for', nrow(pathrows), 
                                            'path/rows'))
}
