#' Generate random sample polygons from a raster layer, optionally with 
#' stratification
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
#' @param strata (optional) a \code{RasterLayer} of integers giving the strata 
#' of each pixel.
#' @param fields a list of fields to include in the output 
#' \code{SpatialPolygonsDataFrame} (such as a "class" field if you will be 
#' digitizing classes).
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
#' training_polys <- sample_observed(L5TSR_1986_b1, 30,
#'                                   side=6*xres(L5TSR_1986_b1))
#' plot(L5TSR_1986_b1)
#' plot(training_polys, add=TRUE)
#' }
sample_observed <- function(x, size, strata=NULL, side=xres(x), fields=c(), 
                            validate=TRUE, validate.layer=1, exp=5) {
    # Don't expand if no validation is being performed
    if (!validate) {exp <- 1}

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
        strat_sample <- sampleStratified(strata, size, exp=exp)
        cell_nums <- strat_sample[, 1]
        strataids <- strat_sample[, 2]
    } else {
        stratified <- FALSE
        cell_nums <- sampleInt(ncell(x), size * exp)
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
        if (stratified) {
            freqs <- table(polys$strata)
            for (n in 1:length(freqs)) {
                if (freqs[n] > size) {
                    this_class <- names(freqs)[n]
                    these_polys <- which(polys$strata == this_class)
                    drop_polys <- these_polys[sampleInt(length(these_polys), freqs[n] - size)]
                    polys <- polys[-drop_polys, ]
                }
            }
        } else {
            if (length(polys) > size) {
                polys <- polys[sampleInt(length(polys), size), ]
            } else if (length(polys) < size) {
                warning('length of polys < size. Try increasing exp.')
            }
        }
        # Ensure IDs are stills sequential
        polys$ID <- seq(1, length(polys))
    }

    return(polys)
}
