#' Performs a rough comparison of two proj4strings to see if they match
#'
#' Compares the proj, ellps, zone (if applicable), units, and datum tags in two 
#' proj4strings to determine if two projections match. Requires proj and ellps 
#' to be present. If present, zone, units, and datum must match in both 
#' strings.
#'
#' @export
#' @importFrom rgdal CRSargs
#' @importFrom stringr str_extract
#' @param x a proj4string to compare with \code{y}
#' @param y a proj4string to compare with \code{x}
#' @return TRUE or FALSE depending on if the projections match
proj4comp <- function(x, y) {
    if (grepl('+init=', x)) {
        x <- CRSargs(CRS(x))
    }
    if (grepl('+init=', y)) {
        y <- CRSargs(CRS(y))
    }

    x_proj <- str_extract(x, '+proj=[a-zA-Z0-9]*')
    y_proj <- str_extract(y, '+proj=[a-zA-Z0-9]*')
    if (sum(is.na(c(x_proj, y_proj))) > 0) {
        stop('proj string is missing')
    } else {
        if (x_proj != y_proj) {
            return(FALSE)
        }
    }

    x_ellps <- str_extract(x, '+ellps=[a-zA-Z0-9]*')
    y_ellps <- str_extract(y, '+ellps=[a-zA-Z0-9]*')
    if (sum(is.na(c(x_ellps, y_ellps))) > 0) {
        stop('ellps string is missing')
    } else if (x_ellps != y_ellps) {
        return(FALSE)
    }


    if (grepl('utm', tolower(x_proj))) {
        x_zone <- str_extract(x, '+zone=[a-zA-Z0-9]*')
        y_zone <- str_extract(y, '+zone=[a-zA-Z0-9]*')
        if (sum(is.na(c(x_zone, y_zone))) > 0) {
            stop('utm zone must be specified for utm projections')
        } else if (x_zone != y_zone) {
            return(FALSE)
        }
        x_south <- str_extract(tolower(x), '+south')
        y_south <- str_extract(tolower(y), '+south')
        if (!is.na(x_south) || !is.na(y_south)) {
            # Get here if one or more of x_south and y_south is not NA
            if (xor(is.na(x_south), is.na(y_south)) || (x_south != y_south)) {
                # Get here if ONLY one of x_south and y_south is NA, or, if 
                # x_south and y_south are both not NA and are both not equal
                return(FALSE)
            }
        }
    }

    x_datum <- str_extract(x, '+datum=[a-zA-Z0-9]*')
    y_datum <- str_extract(y, '+datum=[a-zA-Z0-9]*')
    if ((!is.na(x_datum) && !is.na(y_datum)) & (x_datum != y_datum)) {
        return(FALSE)
    }
    
    x_units <- str_extract(x, '+units=[a-zA-Z0-9]*')
    y_units <- str_extract(y, '+units=[a-zA-Z0-9]*')
    if ((!is.na(x_units) && !is.na(y_units)) & (x_units != y_units)) {
        return(FALSE)
    }

    return(TRUE)
}
