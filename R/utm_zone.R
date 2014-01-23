#' Given a spatial object, calculate the UTM zone of the centroid
#'
#' For a line or polygon, the UTM zone of the centroid is given, after 
#' reprojecting the object into WGS-84.
#'
#' Based on the code on gis.stackexchange.com at http://bit.ly/17SdcuN.
#'
#' @export
#' @param long a longitude (with western hemisphere longitudes negative)
#' @param lat a latitude (with southern hemisphere latitudes negative)
#' @param proj4string if FALSE (default) return the UTM zone as a string (for 
#' example "34S" for UTM Zone 34 South). If TRUE, return a proj4string using 
#' the EPSG code as an initialization string.
#' @examples
#' utm_zone(45, 10)
#' utm_zone(45, -10)
#' utm_zone(45, 10, proj4string=TRUE)
utm_zone <- function(long, lat, proj4string=FALSE) {
    if (long < -180 || long > 180) {
        stop("longitude must be between -180 and 180")
    }
    if (lat < -90 || lat > 90) {
        stop("latitude must be between -90 and 90")
    }

    zone_num <- floor((long + 180)/6) + 1
    if (lat >= 56.0 && lat < 64.0 && long >= 3.0 && long < 12.0) {
        zone_num <- 32
    }

    # Special zone_nums for Svalbard
    if (lat >= 72.0 && lat < 84.0) {
        if (long >= 0.0 && long < 9.0) {
            zone_num <- 31
        } else if (long >= 9.0 && long < 21.0) {
            zone_num <- 33
        } else if (long >= 21.0 && long < 33.0) {
            zone_num <- 35
        } else if (long >= 33.0 && long < 42.0) {
            zone_num <- 37
        }
    }

    if (lat >= 0) {
        ns <- 'N'
    } else {
        ns <- 'S'
    }

    if (proj4string) {
        if (ns == 'N') {
            return(paste0('+init=epsg:326', sprintf('%02i', zone_num)))
        } else {
            return(paste0('+init=epsg:327', sprintf('%02i', zone_num)))
        }
    } else {
        return(paste0(zone_num, ns))
    }
}
