#' Assign class labels to classification file
#'
#' @export
#' @importFrom rgdal writeGDAL
#' @importFrom sp SpatialPixelsDataFrame
#' @param x classified image as \code{RasterLayer}
#' @param cls two column matrix, where the first column is the class codes 
#' (integers) and the second column is the class names
#' @param outfile the filename to use for the output
#' @examples
#' #TODO: Add examples
color_image <- function(x, cls, outfile) {
    # Reclassify image so it is coded from 0 - length(cls[1]). ENVI and 
    # other classification file formats rely on the codes being sequential, 
    # starting at zero.
    x <- reclassify(x, cbind(cls[, 1], seq(0, length(cls[, 1]) - 1)))
    x.sp <- as(x, "SpatialPixelsDataFrame")
    cls_colors <- t(col2rgb(cls[, 1]))
    # Select appropriate data type and missing value tag depending on data 
    # attributes.
    if (max(x.sp$layer) > 254) {
        gdaltype <- 'Int16'
        gdalmvFlag <- -32768
    } else {
        gdaltype <- 'Byte'
        gdalmvFlag <- 255
    }
    writeGDAL(x.sp, outfile, drivername="ENVI", type=gdaltype,
              colorTables=list(cls_colors), catNames=list(as.character(cls[, 2])),
              mvFlag=gdalmvFlag)
}
