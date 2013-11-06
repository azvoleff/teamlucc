#' Topographic correction for satellite imagery
#'
#' Perform topographic correction using a number of alternative methods. This 
#' code is taken directly from the \code{landsat} package by Sarah Goslee. The 
#' code has been altered from the \code{landsat} version to allow the option of 
#' using a sample of pixels for calculation of k in the Minnaert correction 
#' (useful when dealing with large images).
#' 
#' See the help page for \code{topocorr} in the \code{landsat} package for 
#' details on the parameters.
#'
#' @export
#' @import spatial.tools
#' @importFrom landsat slopeasp
#' @param DEM DEM as \code{RasterLayer}
#' @param EWkernel kernel to use for East-West gradient calculation
#' @param NSkernel kernel to use for North-South gradient calculation
#' @param smoothing positive integer for smoothing. 1 means no smoothing.
#' @param filename output file for 2-band slope and aspect layer stack 
#' (optional)
#' @param overwrite whether to overwrite output file if it exists
#' @param sampleindices (optional) row-major indices of sample pixels to use in 
#' regression models used for some topographic correction methods (like 
#' Minnaert). Useful when handling very large images. See
#' \code{\link{gridsample}} for one method of calculating these indices.
#' @return RasterBrick with two layers: 'slope' and 'aspect'
#' @references
#' Sarah Goslee. Analyzing Remote Sensing Data in {R}: The {landsat} Package.  
#' Journal of Statistical Software, 2011, 43:4, pg 1--25.  
#' http://www.jstatsoft.org/v43/i04/
topocorr_samp <- function(x, slope, aspect, sunelev, sunazimuth, method="cosine", 
                          na.value=NA, GRASS.aspect=FALSE, IL.epsilon=0.000001,
                          sampleindices=NULL) {
    stop('topocorr_samp is not yet supported. try using method="minnaert_full")')
    ## aspect may be GRASS output: counterclockwise from east
    ## or nonGRASS output: clockwise from north
    ## require the latter for further calculations
    ## because sunazimuth is always measured clockwise from north
    if(GRASS.aspect) {
        aspect <- as.matrix(aspect)
        aspect <- -1 * aspect + 90
        aspect <- (aspect + 360) %% 360
    }

    # all inputs are in degrees, but we need radians
    slope <- (pi/180) * as.matrix(slope)
    aspect <- (pi/180) * as.matrix(aspect)
    sunzenith <- (pi/180) * (90 - sunelev)
    sunazimuth <- (pi/180) * sunazimuth

    x.orig <- x
    x <- as.matrix(x)
    x[x == na.value] <- NA

    IL <- cos(slope) * cos(sunzenith) + sin(slope) * sin(sunzenith) * cos(sunazimuth - aspect)
    IL[IL == 0] <- IL.epsilon

    METHODS <- c("cosine", "improvedcosine", "minnaert", "minslope", 
                 "ccorrection", "gamma", "SCS", "illumination")
    method <- pmatch(method, METHODS)
    if (is.na(method)) 
        stop("invalid method")
    if (method == -1) 
        stop("ambiguous method")

    if(method == 1){
        ## Cosine method
        xout <- x * (cos(sunzenith)/IL)
    } else if(method == 2) {
        ## Improved cosine method
        ILmean <- mean(as.vector(IL), na.rm=TRUE)
        xout <- x + (x * (ILmean - IL)/ILmean)
    } else if(method == 3) {
        ## Minnaert
        ## K is between 0 and 1
        ## only use points with greater than 5% slope
        targetslope <- atan(.05)

        if(all(x[slope >= targetslope] < 0, na.rm=TRUE)) {
            K <- 1
        } else {
            if (usesample) {
                horizcells <- 10
                vertcells <- 10
                nsamp <- 100000 / (horizcells * vertcells)
                x_samp <- gridsample(x[slope >= targetslope], nsamp=nsamp, 
                                     horizcells=10, vertcells=10, 
                                     returnindices=TRUE)
                K <- data.frame(y = as.vector(x_samp$value),
                                x = as.vector(IL[slope >= 
                                              targetslope][x_samp$index]) / 
                                              cos(sunzenith))
            } else {
                K <- data.frame(y = as.vector(x[slope >= targetslope]), x = 
                                as.vector(IL[slope >= 
                                          targetslope])/cos(sunzenith))
            }
            # IL can be <=0 under certain conditions
            # but that makes it impossible to take log10 so remove those 
            # elements
            K <- K[!apply(K, 1, function(x)any(is.na(x))),]
            K <- K[K$x > 0, ]
            K <- K[K$y > 0, ]

            K <- lm(log10(K$y) ~ log10(K$x))
            K <- coefficients(K)[[2]] # need slope
            if(K > 1) K <- 1
            if(K < 0) K <- 0
        }

        xout <- x * (cos(sunzenith)/IL) ^ K
    } else if(method == 4) {
        ## Minnaert with slope
        ## K is between 0 and 1
        ## only use points with greater than 5% slope
        targetslope <- atan(.05)

        if(all(x[slope >= targetslope] < 0, na.rm=TRUE)) {
            K <- 1
        } else {
            if (usesample) {
                horizcells <- 10
                vertcells <- 10
                nsamp <- 100000 / (horizcells * vertcells)
                x_samp <- gridsample(x[slope >= targetslope], nsamp=nsamp, 
                                     horizcells=10, vertcells=10, 
                                     returnindices=TRUE)
                K <- data.frame(y = as.vector(x_samp$value),
                                x = as.vector(IL[slope >= 
                                              targetslope][x_samp$index]) / 
                                              cos(sunzenith))
            } else {
                K <- data.frame(y=as.vector(x[slope >= targetslope]), 
                                x=as.vector(IL[slope >= 
                                            targetslope])/cos(sunzenith))
            }
            # IL can be <=0 under certain conditions
            # but that makes it impossible to take log10 so remove those elements
            K <- K[!apply(K, 1, function(x) any(is.na(x))),]
            K <- K[K$x > 0, ]
            K <- K[K$y > 0, ]

            K <- lm(log10(K$y) ~ log10(K$x))
            K <- coefficients(K)[[2]] # need slope
            if(K > 1) K <- 1
            if(K < 0) K <- 0
        }
        xout <- x * cos(slope) * (cos(sunzenith) / (IL * cos(slope))) ^ K
    } else if(method == 5) {
        ## C correction
        x_samp <- gridsample(x, nsamp=nsamp, horizcells=10, vertcells=10, 
                             returnindices=TRUE)
        band.lm <- lm(as.vector(x_samp$value) ~ as.vector(IL[x_samp$index]))
        C <- coefficients(band.lm)[[1]]/coefficients(band.lm)[[2]]

        xout <- x * (cos(sunzenith) + C) / (IL + C)
    } else if(method == 6) {
        ## Gamma
        ## assumes zenith viewing angle
        viewterrain <- pi/2 - slope
        xout <- x * (cos(sunzenith) + cos(pi / 2)) / (IL + cos(viewterrain))
    } else if(method == 7) {
        ## SCS method from GZ2009
        xout <- x * (cos(sunzenith) * cos(slope))/IL
    } else if(method == 8) {
        ## illumination only
        xout <- IL
    }

    ## if slope is zero, reflectance does not change
    if(method != 8) 
        xout[slope == 0 & !is.na(slope)] <- x[slope == 0 & !is.na(slope)]

    ## if x was a SpatialGridDataFrame, return an object of the same class
    if(class(x.orig) == "SpatialGridDataFrame") {
        x.orig@data[,1] <- as.vector(xout)
        xout <- x.orig
    }

    return(xout)
}

