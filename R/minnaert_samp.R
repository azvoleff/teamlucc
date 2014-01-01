#' @import mgcv
.calc_k_table <- function(x, IL, slope, sampleindices, slopeclass,
                          coverclass, sunzenith) {
    if (!is.null(sampleindices)) {
        K <- data.frame(x=x[sampleindices],
                        IL=IL[sampleindices], 
                        slope=slope[sampleindices])
        # Remember that the sample indices are row-major (as they were drawn 
        # for a RasterLayer), so the coverclass matrix needs to be transposed 
        # as it is stored in column-major order
        coverclass <- t(coverclass)[sampleindices]
    } else {
        K <- data.frame(x=getValues(x), IL=getValues(IL), 
                        slope=getValues(slope))
    }

    ## K is between 0 and 1
    # IL can be <=0 under certain conditions
    # but that makes it impossible to take log10 so remove those elements
    K <- K[coverclass, ]
    K <- K[!apply(K, 1, function(rowvals) any(is.na(rowvals))),]
    K <- K[K$x > 0, ]
    K <- K[K$IL > 0, ]

    k_table <- data.frame(matrix(0, nrow=length(slopeclass) - 1, ncol=3))
    names(k_table) <- c('midpoint', 'n', 'k')
    k_table$midpoint <- diff(slopeclass)/2 + slopeclass[1:length(slopeclass) - 1]

    # don't use slopes outside slopeclass range
    K.cut <- as.numeric(cut(K$slope, slopeclass))
    if(nrow(k_table) != length(table(K.cut))) {
        stop("slopeclass is inappropriate for these data (empty classes)\n")
    }
    k_table$n <- table(K.cut)

    for(i in sort(unique(K.cut[!is.na(K.cut)]))) {
        k_table$k[i] <- coefficients(lm(log10(K$x)[K.cut == i] ~ 
                                        log10(K$IL/cos(sunzenith))[K.cut == i]))[[2]]
    }

    return(k_table)
}

#' Topographic correction for satellite imagery using Minnaert method
#'
#' Perform topographic correction using the Minnaert method. This code is 
#' modified from the code in the \code{landsat} package written by Sarah 
#' Goslee.  This version of the code has been altered from the \code{landsat} 
#' version to allow the option of using a sample of pixels for calculation of k 
#' in the Minnaert correction (useful when dealing with large images).
#' 
#' See the help page for \code{minnaert} in the \code{landsat} package for 
#' additional details on the parameters.
#'
#' @export
#' @import spatial.tools
#' @param x image as a \code{RasterLayer}
#' @param slope the slope as a \code{RasterLayer}
#' @param aspect the aspect as a \code{RasterLayer}
#' @param sunelev sun elevation in degrees
#' @param sunazimuth sun azimuth in degrees
#' @param na.value the value used to code no data values
#' @param GRASS.aspect is aspect defined according to GRASS convetion
#' @param IL.epsilon a small amount to add to calculated illumination values 
#' that are equal to zero to avoid division by zero resulting in Inf values
#' @param slopeclass the slope classes to calculate k for
#' @param coverclass used to calculate k for specific cover class (optional)
#' @param sampleindices (optional) row-major indices of sample pixels to use in 
#' the calculation of k values for the Minnaert correction. See
#' \code{\link{gridsample}}.
#' @return \code{RasterLayer} with topographically corrected data
#' @author Sarah Goslee and Alex Zvoleff
#' @references
#' Sarah Goslee. Analyzing Remote Sensing Data in {R}: The {landsat} Package.  
#' Journal of Statistical Software, 2011, 43:4, pg 1--25.  
#' http://www.jstatsoft.org/v43/i04/
minnaert_samp <- function(x, slope, aspect, sunelev, sunazimuth,
                          na.value=NA, GRASS.aspect=FALSE, IL.epsilon=0.000001,
                          slopeclass=c(1, 5, 10, 15, 20, 25, 30, 45), 
                          coverclass=NULL, sampleindices=NULL) {
    ## aspect may be GRASS output: counterclockwise from east
    ## or nonGRASS output: clockwise from north
    ## require the latter for further calculations
    ## because sunazimuth is always measured clockwise from north
    if(GRASS.aspect) {
        aspect <- -1 * aspect + 90
        aspect <- (aspect + 360) %% 360
    }

    x[x == na.value] <- NA

    # all inputs are in degrees, but we need radians
    sunzenith <- (pi/180) * (90 - sunelev)
    sunazimuth <- (pi/180) * sunazimuth
    slope <- (pi/180) * slope
    aspect <- (pi/180) * aspect
    slopeclass <- (pi/180) * slopeclass

    IL <- .calc_IL(slope, sunzenith, sunazimuth, aspect, IL.epsilon)
    rm(aspect, sunazimuth)

    if(is.null(coverclass)) 
        coverclass <- matrix(rep(TRUE, length(x)), nrow=nrow(x))

    k_table <- .calc_k_table(x, IL, slope, sampleindices, slopeclass, 
                             coverclass, sunzenith)

    k_model <- gam(k ~ s(midpoint, k=length(k_table$midpoint) - 1), data=k_table)

    names(slope) <- 'midpoint'
    # if slope is greater than modeled range, use maximum of modeled range
    slope[slope > max(slopeclass)] <- max(slopeclass)
    # if slope is less than modeled range, treat it as flat
    slope[slope < min(slopeclass)] <- 0
    print("Made it to prediction line in minnaert_samp")
    K.all <- predict(slope, k_model)
    K.all[K.all > 1] <- 1
    K.all[K.all < 0] <- 0

    # Perform correction
    xout <- x * (cos(sunzenith)/IL) ^ K.all
    # Don't correct flat areas
    xout[K.all == 0 & !is.na(K.all)] <- x[K.all == 0 & !is.na(K.all)]

    list(classcoef=k_table, model=k_model, minnaert=xout, sampleindices=sampleindices)
}
