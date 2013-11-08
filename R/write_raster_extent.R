#' Makes a polygon shapefile of a \code{Raster*} extent.
#'
#' @export
#' @import raster rgdal
#' @param x A \code{Raster*} object
#' @param out_file Filename for the output shapefile
#' @param fields list of fields to include in the output shapefile
write_raster_extent <- function(x, out_file, fields=c('t1_class',
                                                      't2_class')) {
    if ((length(out_file) > 1) || (!grepl('.shp$', out_file, 
                                          ignore.case=TRUE))) {
        stop("out_file must be a file path ending in '.shp'")
    }
    img_ext <- extent(x)
    # Need to create the shapefile with a polygon in it as R doesn't have an 
    # easy way to create and setup an empty shapefile. So make a polygon with 
    # the image extent to use as the initial polygon in the training data 
    # shapefile.
    img_ext_poly <- Polygons(list(Polygon(cbind(x=c(img_ext@xmin, img_ext@xmax, 
                                                    img_ext@xmax, img_ext@xmin, 
                                                    img_ext@xmin),
                                                y=c(img_ext@ymin, img_ext@ymin, 
                                                    img_ext@ymax, img_ext@ymax, 
                                                    img_ext@ymin)))), ID=1)
    Sr <- SpatialPolygons(list(img_ext_poly), 
                          proj4string=CRS(proj4string(x)))
    # Create a data.frame with empty rows for each field in field, and with the 
    # 'Notes' field set to 'Extent'
    out_data <- data.frame(t(c(rep('', length(fields)), 'Extent')))
    names(out_data) <- c(fields, 'Notes')
    # Now create the training data SpatialPolygonsDataFrame instance
    train_df <- SpatialPolygonsDataFrame(Sr, data=data.frame(out_data))
    layer_name <- gsub('.shp$', '', basename(out_file), ignore.case=TRUE)
    writeOGR(train_df, dirname(out_file), layer_name, driver='ESRI Shapefile')
}
