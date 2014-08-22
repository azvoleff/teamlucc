#' Generate random sample polygons from a raster layer, optionally with 
#' stratification
#'
#' Useful for gathering training data for an image classification. With the 
#' default settings, the output polygons will be perfectly aligned with the 
#' pixels in the input raster.
#'
#' @export
#' @importFrom rgdal writeOGR
#' @param x a \code{Raster*}
#' @param size the sample size (number of sample polygons to return)
#' @param side desired length for each side of the sample polygon (units of the 
#' input \code{Raster*}, usually meters)
#' @param strata (optional) a \code{RasterLayer} of integers giving the strata 
#' of each pixel (for example, a classified image)
#' @param fields a list of fields to include in the output 
#' \code{SpatialPolygonsDataFrame} (such as a "class" field if you will be 
#' digitizing classes).
#' @param na.rm whether to remove pixels with NA values from the sample
#' @param exp multiplier used to draw larger initial sample to account for the 
#' loss of sample polygons lost because they contain NAs, and, for stratified 
#' sampling, to account for classes that occur very infrequently in the data.  
#' Increase this value if the final sample has fewer sample polygons than 
#' desired.
#' @return a \code{SpatialPolygonsDataFrame}
#' @examples
#' \dontrun{
#' set.seed(0)
#' L5TSR_1986_b1 <- raster(L5TSR_1986, layer=1)
#' training_polys <- sample_raster(L5TSR_1986_b1, 30,
#'                                   side=6*xres(L5TSR_1986_b1))
#' plot(L5TSR_1986_b1)
#' plot(training_polys, add=TRUE)
#' }
sample_raster <- function(x, size, strata=NULL, side=xres(x), fields=c(), 
                          na.rm=TRUE, exp=5) {
    if (!is.null(strata)) {
        stratified <- TRUE
        if (proj4string(strata) != proj4string(x)) {
            stop('x and strata must have the same coordinate system')
        }
        if (!identical(extent(strata), extent(x))) {
            stop('x and strata must have the same extent')
        }
        if (!identical(res(strata), res(x))) {
            stop('x and strata must have the same resolution')
        }
        strat_sample <- sampleStratified(strata, size, exp=exp, na.rm=na.rm)
        cell_nums <- strat_sample[, 1]
        strataids <- strat_sample[, 2]
    } else {
        stratified <- FALSE
        cell_nums <- sampleRandom(x, size, exp=exp, na.rm=na.rm)
    }

    xy <- xyFromCell(x, cell_nums)
    # Convert from cell-center coordinates to ul corner of cell coordinates
    xy[, 1] <- xy[, 1] - xres(x)/2
    xy[, 2] <- xy[, 2] + yres(x)/2
    # Coordinate order is: ll, lr, ur, ul, ll. Need to end with ll to close the 
    # polygon
    xcoords <- cbind(xy[, 1],
                     xy[, 1] + side,
                     xy[, 1] + side,
                     xy[, 1],
                     xy[, 1])
    ycoords <- cbind(xy[, 2] - side,
                     xy[, 2] - side,
                     xy[, 2],
                     xy[, 2],
                     xy[, 2] - side)
    xycoords <- array(cbind(xcoords, ycoords),
                      dim=c(nrow(xcoords), ncol(xcoords), 2))

    # TODO: Add check to ensure no polygons overlap, and warn if they fall 
    # outside raster extent

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
    if (stratified) {
        out_data <- data.frame(ID=names(Sr), strata=strataids)
    } else {
        out_data <- data.frame(ID=names(Sr))
    }
    for (field in fields) {
        out_data <- cbind(out_data, rep('', nrow(out_data)))
        names(out_data)[ncol(out_data)] <- field
    }
    # Finally, convert to a SpatialPolygonsDataFrame
    polys <- SpatialPolygonsDataFrame(Sr, data=out_data)

    if (nrow(polys) < size) {
        warning('length of polys < size. Try increasing exp.')
    }

    return(polys)
}
