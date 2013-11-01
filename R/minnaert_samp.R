#' Topographic correction for satellite imagery
#'
#' Perform topographic correction using the Minnaert method. This code is taken 
#' directly from the \code{landsat} package by Sarah Goslee.  The code has been 
#' altered from the \code{landsat} version to allow the option of using a 
#' sample of pixels for calculation of k in the Minnaert correction (useful 
#' when dealing with large images).
#' 
#' See the help page for \code{minnaert} in the \code{landsat} package for 
#' details on the parameters.
#'
#' @export
#' @import spatial.tools
#' @import mgcv
#' @importFrom landsat slopeasp
#' @param DEM DEM as \code{RasterLayer}
#' @param EWkernel kernel to use for East-West gradient calculation
#' @param NSkernel kernel to use for North-South gradient calculation
#' @param smoothing positive integer for smoothing. 1 means no smoothing.
#' @param filename output file for 2-band slope and aspect layer stack 
#' (optional)
#' @param overwrite whether to overwrite output file if it exists
#' @return RasterBrick with two layers: 'slope' and 'aspect'
#' @references
#' Sarah Goslee. Analyzing Remote Sensing Data in {R}: The {landsat} Package.  
#' Journal of Statistical Software, 2011, 43:4, pg 1--25.  
#' http://www.jstatsoft.org/v43/i04/
minnaert_samp <- function(x, slope, aspect, sunelev, sunazimuth, na.value=NA, 
                     GRASS.aspect=FALSE, IL.epsilon=0.000001,
                     slopeclass = c(1, 5, 10, 15, 20, 25, 30, 45), 
                     usesample=FALSE, coverclass) {
    # IL.epsilon: if IL == 0, the corrected value is Inf (division by zero)
    # adding a tiny increment eliminates the Inf

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
    sloper <- (pi/180) * as.matrix(slope)
    sloped <- as.matrix(slope)
    aspect <- (pi/180) * as.matrix(aspect)
    sunzenith <- (pi/180) * (90 - sunelev)
    sunazimuth <- (pi/180) * sunazimuth

    x.orig <- x
    x <- as.matrix(x)
    x[x == na.value] <- NA

    IL <- cos(sloper) * cos(sunzenith) + sin(sloper) * sin(sunzenith) * cos(sunazimuth - aspect)
    IL[IL == 0] <- IL.epsilon

    if(missing(coverclass)) 
        coverclass <- rep(TRUE, length(as.vector(x)))

    ## Minnaert
    ## K is between 0 and 1
    # IL can be <=0 under certain conditions
    # but that makes it impossible to take log10 so remove those elements
    if (usesample) {
        horizcells <- 10
        vertcells <- 10
        nsamp <- 100000 / (horizcells * vertcells)
        x_samp <- gridsample(x, nsamp=nsamp, horizcells=10, vertcells=10, 
                             returnindices=TRUE)
        K <- data.frame(x = as.vector(x_samp$value), IL = 
                        as.vector(IL[x_samp$index]), 
                        slope=as.vector(sloped[x_samp$index]))
    } else {
        K <- data.frame(x = as.vector(x), IL = as.vector(IL), 
                        slope=as.vector(sloped))
    }
    K <- K[coverclass, ]
    K <- K[!apply(K, 1, function(x)any(is.na(x))),]
    K <- K[K$x > 0, ]
    K <- K[K$IL > 0, ]

    ## calculate overall K value; only use points with greater than 5% slope
    targetslope <- (180/pi) * atan(.05)
    allcoef <- coefficients(lm(log10(K$x)[K$slope >= targetslope] ~ log10(K$IL/cos(sunzenith))[K$slope >= targetslope]))[[2]]

    results <- data.frame(matrix(0, nrow=length(slopeclass)-1, ncol=3))
    colnames(results) <- c("midpoint", "n", "k")
    results[,1] <- diff(slopeclass)/2 + slopeclass[1:length(slopeclass)-1]

    K.cut <- as.numeric(cut(K$slope, slopeclass)) # don't use slopes outside slopeclass range
    if(nrow(results) != length(table(K.cut))) {
        stop("slopeclass is inappropriate for these data (empty classes)\n")
    }
    results[,2] <- table(K.cut)

    #            sapply(unique(K.cut), function(i)coefficients(lm(log10(K$x)[K.cut == i] ~ log10(K$IL/cos(sunzenith))[K.cut == i]))[[2]])
    for(i in sort(unique(K.cut[!is.na(K.cut)]))) {
        results[i, 3] <- coefficients(lm(log10(K$x)[K.cut == i] ~ 
                                         log10(K$IL/cos(sunzenith))[K.cut == 
                                                                    i]))[[2]]
    }

    model <- with(results, gam(k ~ s(midpoint, k=length(midpoint)-1)))

    K.all <- data.frame(midpoint = as.vector(as.matrix(slope)))
    K.all[K.all > max(slopeclass)] <- max(slopeclass) # if slope is greater than modeled range, use maximum of modeled range
    K.all[K.all < min(slopeclass)] <- 0 # if slope is less than modeled range, treat it as flat
    K.all <- predict(model, newdata=K.all)
    K.all[K.all > 1] <- 1
    K.all[K.all < 0] <- 0

    xout <- as.vector(as.matrix(x)) * (cos(sunzenith)/as.vector(as.matrix(IL))) ^ K.all
    xout[K.all == 0 & !is.na(K.all)] <- as.vector(as.matrix(x))[K.all == 0 & !is.na(K.all)] # don't correct flat areas

    ## if x was a SpatialGridDataFrame, return an object of the same class
    if(class(x.orig) == "SpatialGridDataFrame") {
        x.orig@data[,1] <- as.vector(xout)
        xout <- x.orig
    }

    list(allcoef=allcoef, classcoef=results, model=model, minnaert=xout)

}

