#' Generate random sample polygons from a raster layer
#'
#' Useful for gathering training data for an image classification. With the 
#' default settings, the output polygons will be perfectly aligned with the 
#' pixels in the input raster.
#'
#' The \code{validate} parameter allows excluding no data areas from the output 
#' polygons (though the simple method this function uses to enforce this means 
#' that pixels within the border areas of the image will have a slightly lower 
#' probability of selection than those in the interior).
#'
#' @export
#' @importFrom rgdal writeOGR
#' @param x a \code{Raster*}
#' @param size the sample size (number of sample polygons to return)
#' @param side the length each side of the sample polygon
#' @param fields a list of fields to include in the output 
#' \code{SpatialPolygonsDataFrame}
#' @param out_file (optional) shapefile to save to save output
#' @param validate whether to check that all sample polygons lie within image 
#' area (defined as having no NAs within the polygon)
#' @param validate.layer the index within \code{x} of the layer that should be 
#' used for validation
#' @param exp multiplier used to draw larger initial sample to account for the 
#' loss of sample polygons during validation. Increase this value if the final 
#' sample have fewer sample polygons than desired \code{size}
#' @return a \code{SpatialPolygonsDataFrame} with \code{n} polygons
#' @examples
#' \dontrun{
#' set.seed(0)
#' L5TSR_1986_b1 <- raster(L5TSR_1986, layer=1)
#' training_polys <- training_sample(L5TSR_1986_b1, 30,
#'                                   side=6*xres(L5TSR_1986_b1))
#' plot(L5TSR_1986_b1)
#' plot(training_polys, add=TRUE)
#' }
training_sample <- function(x, size, side=2*xres(x), fields=c(), out_file=NULL, 
                          validate=TRUE, validate.layer=1, exp=2) {
    # Don't expand if no validation is being performed
    if (!validate) {exp <- 1}

    cell_nums <- sampleInt(ncell(x), size * exp)

    xy <- xyFromCell(x, cell_nums)
    # Convert from cell-center coordinates to ul corner of cell coordinates
    xy[, 1] <- xy[, 1] - xres(x)/2
    xy[, 2] <- xy[, 2] + yres(x)/2
    # Coordinate order is: ll, lr, ur, ul, ll. Need to end with ll to close the 
    # polygon
    xcoords <- cbind(xy[, 1] - side/2,
                     xy[, 1] + side/2,
                     xy[, 1] + side/2,
                     xy[, 1] - side/2,
                     xy[, 1] - side/2)
    ycoords <- cbind(xy[, 2] - side/2,
                     xy[, 2] - side/2,
                     xy[, 2] + side/2,
                     xy[, 2] + side/2,
                     xy[, 2] - side/2)
    xycoords <- array(cbind(xcoords, ycoords),
                      dim=c(nrow(xcoords), ncol(xcoords), 2))

    # TODO: Add check to ensure no polygons overlap, and warn if they fall 
    # outside rater extent

    # Function to make individual polygons for each sample area
    make_Polygon <- function(slice) {
        list(Polygon(cbind(x=slice[, 1], y=slice[, 2])))
    }
    polys <- apply(xycoords, c(1), make_Polygon)
    # Now convert the list of Polygon objects to a list of Polygons objects 
    # (notice the trailing "s" in "Polygons")
    polys <- mapply(function(poly, ID) Polygons(poly, ID=ID),
                    polys, seq(1, length(polys)))
    Sr <- SpatialPolygons(polys, proj4string=CRS(proj4string(x)))
    out_data <- data.frame(ID=names(Sr))
    for (field in fields) {
        out_data <- cbind(out_data, rep('', nrow(out_data)))
        names(out_data)[ncol(out_data)] <- field
    }
    # Finally, convert to a SpatialPolygonsDataFrame
    polys <- SpatialPolygonsDataFrame(Sr, data=out_data)

    if (validate) {
        if (nlayers(x) > 1) {
            have_NAs <- extract(raster(x, layer=validate.layer), polys,
                                fun=function(windo, ...) any(is.na(windo)))
        } else {
            have_NAs <- extract(x, polys,
                                fun=function(windo, ...) any(is.na(windo)))
        }
    
        # Exclude polys with NAs (as they fall in no data areas)
        polys <- polys[!(have_NAs), ]

        # If there are too still too many polys remaining in the sample, cut sample 
        # down to desired size,
        if (length(polys) > size) {
            polys <- polys[sampleInt(length(polys), size), ]
        } else if (length(polys) < size) {
            warning('length of polys < size. Try increasing exp.')
        }
        # Ensure IDs are stills sequential
        polys$ID <- seq(1, length(polys))
    }

    if (!is.null(out_file)) {
        layer_name <- gsub('.shp$', '', basename(out_file), ignore.case=TRUE)
        writeOGR(polys, dirname(out_file), layer_name, driver='ESRI Shapefile')
    }
    
    return(polys)
}
