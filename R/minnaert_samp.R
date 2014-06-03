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
        if(!is.null(coverclass)) coverclass <- t(coverclass)[sampleindices]
    } else {
        K <- data.frame(x=getValues(x), IL=getValues(IL), 
                        slope=getValues(slope))
    }

    if(!is.null(coverclass)) {
        K <- K[coverclass, ]
    }

    ## K is between 0 and 1
    # IL can be <=0 under certain conditions
    # but that makes it impossible to take log10 so remove those elements
    K <- K[!apply(K, 1, function(rowvals) any(is.na(rowvals))), ]
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

# Function to combine classes defined by the upper limits 'lims' with their 
# smallest neighbors until each class has at least n members. 'counts' stores 
# the number of members in each class.
clean_intervals <- function(counts, lims, n) {
    while(min(counts) < n) {
        # The "rev" below is so that the classes at the end are combined first 
        # (as classes with higher slopes are more more likely to be more rare 
        # and therefore have fewer members)
        min_index <- length(counts) - match(TRUE,
                                            rev(counts == min(counts))) + 1
        if (min_index == length(counts)) {
            counts[min_index - 1] <- counts[min_index - 1] + counts[min_index]
            lims <- lims[-(min_index - 1)]
        } else if (min_index == 1) {
            counts[min_index + 1] <- counts[min_index + 1] + counts[min_index]
            lims <- lims[-min_index]
        } else if (counts[min_index - 1] < counts[min_index + 1]) {
            counts[min_index - 1] <- counts[min_index - 1] + counts[min_index]
            lims <- lims[-(min_index - 1)]
        } else {
            # We know counts[min_index - 1] > counts[min_index + 1]
            counts[min_index + 1] <- counts[min_index + 1] + counts[min_index]
            lims <- lims[-min_index]
        }
        counts <- counts[-min_index]
    }
    return(lims)
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
#' @param x image as a \code{RasterLayer}
#' @param slope the slope in radians as a \code{RasterLayer}
#' @param aspect the aspect in radians as a \code{RasterLayer}
#' @param sunelev sun elevation in degrees
#' @param sunazimuth sun azimuth in degrees
#' @param IL.epsilon a small amount to add to calculated illumination values 
#' that are equal to zero to avoid division by zero resulting in Inf values
#' @param slopeclass the slope classes to calculate k for (in radians), or 
#' NULL, in which case an algorithm will be used to choose reasonable defaults 
#' for the given image. If provided, \code{slopeclass} should be a list of 
#' slope class limits. For example: c(1, 5, 10, 15, 20, 25, 30, 45) * (pi/180)
#' @param coverclass used to calculate k for specific cover class (optional)
#' as \code{RasterLayer}
#' @param sampleindices (optional) row-major indices of sample pixels to use in 
#' the calculation of k values for the Minnaert correction. See
#' \code{\link{gridsample}}.
#' @param DN_min minimum allowable pixel value after correction (values less 
#' than \code{DN_min} are set to NA)
#' @param DN_max maximum allowable pixel value after correction (values less 
#' than \code{DN_max} are set to NA)
#' @return \code{RasterLayer} with topographically corrected data
#' @author Sarah Goslee and Alex Zvoleff
#' @references
#' Sarah Goslee. Analyzing Remote Sensing Data in {R}: The {landsat} Package.  
#' Journal of Statistical Software, 2011, 43:4, pg 1--25.  
#' http://www.jstatsoft.org/v43/i04/
minnaert_samp <- function(x, slope, aspect, sunelev, sunazimuth,
                          IL.epsilon=0.000001, slopeclass=NULL, 
                          coverclass=NULL, sampleindices=NULL, DN_min=NULL, 
                          DN_max=NULL) {

    if (is.null(slopeclass)) {
        slopeclass <- c(1, 2, 3, 4, 5, 6, 8, 10, 12,
                        15, 20, 25, 30, 45, 75) * (pi/180)
    }

    if (is.null(sampleindices)) {
        counts <- raster::freq(raster::cut(slope, slopeclass,
                                           include.lowest=TRUE), 
                               useNA='no')
        # Eliminate empty bins:
        slopeclass <- slopeclass[c(TRUE, 1:(length(slopeclass) - 1) %in% counts[, 1])]
        counts <- counts[, 2]
    } else {
        counts <- as.numeric(table(cut(slope[sampleindices], slopeclass, 
                                       include.lowest=TRUE), useNA='no'))
    }
    # The [-1] below is because clean_intervals only needs the upper limits
    slopeclass <- clean_intervals(counts, slopeclass[-1], 100)
    if (length(slopeclass) <= 5) {
        stop('insufficient sample size to develop k model - try changing slopeclass or sampleindices')
    }
    slopeclass <- c(1*pi/180, slopeclass)

    stopifnot(all((slopeclass >= 0) & slopeclass <= pi/2))

    # some inputs are in degrees, but we need radians
    stopifnot((sunelev >= 0) & (sunelev <= 90))
    stopifnot((sunazimuth >= 0) & (sunazimuth <= 360))
    sunzenith <- (pi/180) * (90 - sunelev)
    sunazimuth <- (pi/180) * sunazimuth

    IL <- .calc_IL(slope, aspect, sunzenith, sunazimuth, IL.epsilon)
    rm(aspect, sunazimuth)

    k_table <- .calc_k_table(x, IL, slope, sampleindices, slopeclass, 
                             coverclass, sunzenith)
    
    k_model <- with(k_table, bam(k ~ s(midpoint, k=length(midpoint) - 1), data=k_table))

    # If slope is greater than modeled range, use maximum of modeled range. If 
    # slope is less than modeled range, treat it as flat.
    slopeclass_max <- max(slopeclass)
    slopeclass_min <- min(slopeclass)
    slope <- calc(slope,
                  fun=function(vals) {
                      vals[vals > slopeclass_max] <- slopeclass_max
                      vals[vals < slopeclass_min] <- 0
                      return(vals)
                  })

    names(slope) <- 'midpoint'
    K.all <- predict(slope, k_model)
    K.all <- calc(K.all,
                  fun=function(vals) {
                      vals[vals > 1] <- 1
                      vals[vals < 0] <- 0
                      return(vals)
                  })

    # Perform correction
    xout <- overlay(x, IL, K.all,
                    fun=function(x_vals, IL_vals, K.all_vals) {
                        x_vals * (cos(sunzenith)/IL_vals) ^ K.all_vals
                    })
    # Don't correct flat areas
    xout[K.all == 0 & !is.na(K.all)] <- x[K.all == 0 & !is.na(K.all)]

    if ((!is.null(DN_min)) || (!is.null(DN_max))) {
        xout <- calc(xout, fun=function(vals) {
                        if (!is.null(DN_min)) vals[vals < DN_min] <- NA
                        if (!is.null(DN_max)) vals[vals > DN_max] <- NA
                        return(vals)
                     })
    }

    list(classcoef=k_table, model=k_model, minnaert=xout, sampleindices=sampleindices)
}
