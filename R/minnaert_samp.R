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
#' @import spatial.tools sp mgcv
#' @param x image as a \code{RasterLayer}
#' @param slope the slope as a \code{RasterLayer}
#' @param aspect the aspect as a \code{RasterLayer}
#' @param sunelev sun elevation in degrees
#' @param sunazimuth sun azimuth in degrees
#' @param IL.epsilon a small amount to add to calculated illumination values 
#' that are equal to zero to avoid division by zero resulting in Inf values
#' @param slopeclass the slope classes to calculate k for
#' @param coverclass used to calculate k for specific cover classes (optional)
#' @param usesample whether to use a sample of the data for calculation of k 
#' for Minnaert correction. See \code{\link{gridsample}}.
#' @return \code{RasterLayer} with topographically corrected data
#' @references
#' Sarah Goslee. Analyzing Remote Sensing Data in {R}: The {landsat} Package.  
#' Journal of Statistical Software, 2011, 43:4, pg 1--25.  
#' http://www.jstatsoft.org/v43/i04/
minnaert_samp <- function(x, slope, aspect, sunelev, sunazimuth,
                          IL.epsilon=0.000001, slopeclass=c(1, 5, 10, 15, 20, 
                                                            25, 30, 45), 
                          coverclass=NULL, usesample=FALSE, ...) {
    # all inputs are in degrees, but we need radians
    sunzenith <- (pi/180) * (90 - sunelev)
    sunazimuth <- (pi/180) * sunazimuth
    slope <- (pi/180) * slope
    aspect <- (pi/180) * aspect
    slopeclass <- (pi/180) * slopeclass

    IL <- cos(slope) * cos(sunzenith) + sin(slope) * sin(sunzenith) *
            cos(sunazimuth - aspect)
    IL[IL == 0] <- IL.epsilon

    if(is.null(coverclass)) 
        coverclass <- rep(TRUE, length(x))

    if (usesample) {
        horizcells <- 10
        vertcells <- 10
        nsamp <- 100000 / (horizcells * vertcells)
        x_samp <- gridsample(x, nsamp=nsamp, horizcells=10, vertcells=10, 
                             returnindices=TRUE)
        K <- data.frame(x=as.vector(x_samp$value),
                        IL=getValues(IL)[x_samp$index], 
                        slope=getValues(slope)[x_samp$index])
        coverclass <- coverclass[x_samp$index]
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

    ## calculate overall K value; only use points with greater than 5% slope
    targetslope <- atan(.05)
    allcoef <- coefficients(lm(log10(K$x)[K$slope >= targetslope] ~ log10(K$IL/cos(sunzenith))[K$slope >= targetslope]))[[2]]

    results <- data.frame(matrix(0, nrow=length(slopeclass)-1, ncol=3))
    colnames(results) <- c("midpoint", "n", "k")
    results[,1] <- diff(slopeclass)/2 + slopeclass[1:length(slopeclass)-1]

    # don't use slopes outside slopeclass range
    K.cut <- as.numeric(cut(K$slope, slopeclass))
    if(nrow(results) != length(table(K.cut))) {
        stop("slopeclass is inappropriate for these data (empty classes)\n")
    }
    results[,2] <- table(K.cut)

    for(i in sort(unique(K.cut[!is.na(K.cut)]))) {
        results[i, 3] <- coefficients(lm(log10(K$x)[K.cut == i] ~ 
                                         log10(K$IL/cos(sunzenith))[K.cut == 
                                                                    i]))[[2]]
    }

    model <- with(results, gam(k ~ s(midpoint, k=length(midpoint) - 1)))

    K.all <- data.frame(midpoint=getValues(slope))
    # if slope is greater than modeled range, use maximum of modeled range
    K.all[K.all > max(slopeclass)] <- max(slopeclass)
    # if slope is less than modeled range, treat it as flat
    K.all[K.all < min(slopeclass)] <- 0
    K.all <- predict(model, newdata=K.all)
    K.all[K.all > 1] <- 1
    K.all[K.all < 0] <- 0
    K.all <- raster(matrix(K.all, nrow=nrow(x), byrow=TRUE), template=x)

    # Perform correction
    xout <- x * (cos(sunzenith)/IL) ^ K.all
    xout[K.all == 0 & !is.na(K.all)] <- x[K.all == 0 & !is.na(K.all)] # don't correct flat areas

    list(allcoef=allcoef, classcoef=results, model=model, minnaert=xout)
}

