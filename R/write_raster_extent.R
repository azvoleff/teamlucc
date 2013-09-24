#' Makes a polygon shapefile of image extent.
#'
#' @export
#' @param x A \code{Raster*} object
#' @param out_file Filename for the output shapefile.
write_raster_extent <- function(x, out_file) {
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
    # Now create the training data SpatialPolygonsDataFrame instance
    train_df <- SpatialPolygonsDataFrame(Sr, data=data.frame(Poly_Type='Extent', 
                                                             Notes=''))
    layer_name <- gsub('.shp$', '', basename(out_file), ignore.case=TRUE)
    writeOGR(train_df, dirname(out_file), layer_name, driver='ESRI Shapefile')
}
