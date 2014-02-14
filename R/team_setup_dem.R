#' Setup the DEM mosaic for a particular TEAM site
#'
#' @export
#' @importFrom wrspathrow pathrow_num
#' @importFrom rgdal readOGR writeOGR
#' @importFrom sp spTransform
#' @importFrom rgeos gBuffer gIntersects gUnaryUnion gIntersection
#' @importFrom tools file_path_sans_ext
#' @param dem_path a list of digital elevation models (DEMs) that (when 
#' mosaiced) covers the full extent of all the images in the image_list.
#' @param aoi_file area of interest (AOI) shapefile to use as as bounding box 
#' when selecting DEMs, in a file format readable by \code{readOGR}
#' @param output_path the path to use for the output
#' @param n_cpus the number of CPUs to use for processes that can run in 
#' parallel
#' @param overwrite whether to overwrite existing files (otherwise an error 
#' will be raised)
#' @param crop_to_aoi whether to crop the dem to the supplied AOI, or to the
#' Landsat path/row polygon for that particular path/row
#' @param notify notifier to use (defaults to \code{print} function). See the 
#' \code{notifyR} package for one way of sending notifications from R. The 
#' \code{notify} function should accept a string as the only argument.
#' @examples
#' \dontrun{
#' dem_path <- 'H:/Data/TEAM/VB/Rasters/DEM/ASTER'
#' team_setup_dem(dem_path, "VB", 'H:/Data/TEAM/VB/LCLUC_Analysis/', 
#'                list(c(15,53)))
#' }
team_setup_dem <- function(dem_path, output_path, aoi_file, n_cpus=1, 
                           overwrite=FALSE, crop_to_aoi=FALSE, notify=print) {
    if (!file_test("-d", dem_path)) {
        stop(paste(dem_path, "does not exist"))
    }
    if (!file_test("-d", output_path)) {
        stop(paste(output_path, "does not exist"))
    }

    timer <- Track_time(notify)

    timer <- start_timer(timer, label='Setting up DEMs')

    if (n_cpus > 1) beginCluster(n_cpus)

    # Open datafile with the extents of each CGIAR DEM
    load(file.path(dem_path, 'cgiar_srtm_extents.RData'))

    aoi <- readOGR(dirname(aoi_file), basename(file_path_sans_ext(aoi_file)))
    pathrows <- pathrow_num(aoi, wrs_type=2, wrs_mode='D', as_polys=TRUE)

    writeOGR(pathrows, output_path, 
             paste0(basename(file_path_sans_ext(aoi_file)), '_pathrows'), 
             driver='ESRI Shapefile', overwrite_layer=overwrite)

    png(file.path(output_path, 
        paste0(basename(file_path_sans_ext(aoi_file)), '_pathrows.png')),
        width=900, height=900)
    plot(pathrows, lwd=2)
    aoi_prproj <- spTransform(aoi, CRS(proj4string(pathrows)))
    plot(aoi_prproj, add=TRUE, lty=2, col="#00ff0050", lwd=2)
    text(coordinates(pathrows), labels=paste(pathrows$PATH, pathrows$ROW, 
                                             sep=', '), cex=2)
    dev.off()

    if (crop_to_aoi) {
        pathrows_cropped <- gIntersection(pathrows, aoi_prproj, byid=TRUE)
        row.names(pathrows_cropped) <- row.names(pathrows)
        pathrows <- SpatialPolygonsDataFrame(pathrows_cropped, 
                                             data=pathrows@data)
    }

    # Add a 500 m buffer in UTM coordinate system, as 1) slope calculation 
    # requires a window of pixels, and 2) this buffer also helps avoid missing 
    # pixels on the sides of the image due to slight misalignments from the 
    # reprojection that will occur later. After buffering transform back to 
    # WGS84 to use for preliminary cropping of the dem mosaic. 
    pathrows_utm <- spTransform(pathrows,
                                CRS(utm_zone(pathrows, proj4string=TRUE)))
    pathrows_buffered <- spTransform(gBuffer(pathrows_utm, width=500, byid=TRUE), 
                                 CRS(proj4string(cgiar_srtm_extents)))
    intersecting <- as.logical(gIntersects(cgiar_srtm_extents, 
                                           gUnaryUnion(pathrows_buffered), byid=TRUE))
    if (sum(intersecting) == 0) {
        stop('no intersecting dem extents found')
    } else {
        cgiar_srtm_extents <- cgiar_srtm_extents[intersecting, ]
    }

    dem_list <- file.path(dem_path, cgiar_srtm_extents$filename)
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

    timer <- start_timer(timer, label=paste('Processing', nrow(pathrows),  
                                            'path/rows'))
    for (n in 1:length(pathrows)) {
        pathrow <- pathrows[n, ]
        pathrow_buffered <- pathrows_buffered[n, ]
        pathrow_label <- paste(sprintf('%03i', pathrow@data$PATH), 
                               sprintf('%03i', pathrow@data$ROW), sep='-')

        pathrow_buffered_demproj <- spTransform(pathrow_buffered, 
                                                CRS(proj4string(dem_mosaic)))
        dem_mosaic_crop <- crop(dem_mosaic, pathrow_buffered_demproj)

        timer <- start_timer(timer, label=paste('Reprojecting DEM mosaic crop for', 
                                                pathrow_label))
        # The below lines construct to_ext as the extent the image will be 
        # projected to. This extent must cover the same area as the 
        # dem_mosaic_crop, but must have the same resolution, CRS and origin as 
        # the pathrow
        to_ext <- projectExtent(dem_mosaic_crop,
                                utm_zone(pathrow, proj4string=TRUE))
        to_res <- c(30, 30)
        xmin(to_ext) <- floor(xmin(to_ext) / to_res[1]) * to_res[1]
        xmax(to_ext) <- ceiling(xmax(to_ext) / to_res[1]) * to_res[1]
        ymin(to_ext) <- floor(ymin(to_ext) / to_res[2]) * to_res[2]
        ymax(to_ext) <- ceiling(ymax(to_ext) / to_res[2]) * to_res[2]
        # The below two lines shift the origin to 15, 15 for consistency with 
        # the LEDAPS CDR images.
        xmin(to_ext) <- xmin(to_ext) + to_res[1]/2
        ymax(to_ext) <- ymax(to_ext) + to_res[2]/2
        res(to_ext) <- to_res
        dem_mosaic_crop <- projectRaster(from=dem_mosaic_crop, to=to_ext,
                                         method='bilinear')
        # Crop after the reprojecting to ensure the full crop area is filled.
        if (crop_to_aoi) {
            aoi_demproj <- spTransform(aoi, CRS(proj4string(dem_mosaic_crop)))
            dem_mosaic_crop <- crop(dem_mosaic_crop, aoi_demproj)
        } else {
            pathrow_demproj <- spTransform(pathrow, 
                                           CRS(proj4string(dem_mosaic_crop)))
            dem_mosaic_crop <- crop(dem_mosaic_crop, pathrow_demproj)
        }

        dem_mosaic_filename <- file.path(output_path,
                                         paste0('dem_', pathrow_label, 
                                                '.envi'))
        dem_mosaic_crop <- round(dem_mosaic_crop)
        writeRaster(dem_mosaic_crop, filename=dem_mosaic_filename, 
                    overwrite=overwrite, datatype='INT2S')
        timer <- stop_timer(timer, label=paste('Reprojecting DEM mosaic crop for', 
                                               pathrow_label))

        timer <- start_timer(timer, label=paste('Calculating slope/aspect for', 
                                                pathrow_label))
        slopeaspect_filename <- file.path(output_path,
                                          paste0('slopeaspect_',
                                                 pathrow_label, '.envi'))
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
        timer <- stop_timer(timer, label=paste('Calculating slope/aspect for', 
                                                pathrow_label))
    }
    timer <- stop_timer(timer, label=paste('Processing', nrow(pathrows),  
                                           'path/rows'))

    if (n_cpus > 1) endCluster()

    timer <- stop_timer(timer, label='Setting up DEMs')
}
