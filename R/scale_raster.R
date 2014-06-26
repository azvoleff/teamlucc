#' Scales a \code{Raster*} by a power of a given integer and rounds to nearest 
#' integer
#'
#' Useful for scaling and (optionally) rounding a \code{RasterLayer} to integer 
#' so that a layer can be saved as an integer datatype such as "INT1U", 
#' "INT1S", "INT2" or "INT2S".
#'
#' This function will run in parallel if a parallel backend is registered with 
#' \code{\link{foreach}}.
#'
#' @export scale_raster
#' @import methods
#' @seealso \code{\link{dataType}}
#' @param x a \code{Raster*} object
#' @param power_of raster will be scaled using the highest possible power of 
#' this number
#' @param max_out the scaling factors will be chosen for each layer to ensure 
#' that the maximum and minimum (if minimum is negative) values of each layer 
#' do not exceed \code{max_out}
#' @param do_scaling perform the scaling and return a \code{Raster*} (if 
#' \code{do_scaling} is TRUE) or return a list of scale factors (if 
#' \code{do_scaling} is FALSE)
#' @param round_output whether to round the output to the nearest integer
#' @return a \code{Raster*} if \code{do_scaling} is TRUE, or a list of scaling 
#' factors if \code{do_scaling} is false.
setGeneric("scale_raster", function(x, power_of=10, max_out=32767, 
                                    round_output=TRUE, do_scaling=TRUE) {
    standardGeneric("scale_raster")
})

scale_layer <- function(x, power_of, max_out, round_output, do_scaling) {
    if (!x@data@haveminmax) {
        warning('no stored minimum and maximum values - running setMinMax')
        x <- setMinMax(x)
    }
    layer_max <- max(abs(c(minValue(x), maxValue(x))))
    scale_factor <- power_of ^ floor(log(max_out / layer_max, base=power_of))
    if (do_scaling) {
        x <- calc(x, function(vals, ...) {
                  vals <- vals * scale_factor
                  if (round_output) vals <- round(vals)
                  vals
                  })
        return(x)
    } else {
        return(scale_factor)
    }
}

#' @rdname scale_raster
#' @aliases scale_raster,RasterLayer,ANY-method
setMethod("scale_raster", signature(x="RasterLayer"),
    function(x, power_of, max_out, round_output, do_scaling) {
        ret <- scale_layer(x, power_of, max_out, round_output, do_scaling)
        names(ret) <- names(x)
        return(ret)
    }
)

#' @import foreach
scale_stack_or_brick <- function(x, power_of, max_out, round_output, do_scaling) {
    unscaled_layer=NULL
    if (do_scaling) {
        scale_outputs <- foreach(unscaled_layer=unstack(x), 
                                 .combine='addLayer', .multicombine=TRUE, 
                                 .init=raster(), .packages=c('teamlucc'),
                                 .export=c('scale_layer')) %dopar% {
            scale_output <- scale_layer(unscaled_layer, power_of, max_out, 
                                        round_output, do_scaling)
        }
    } else {
        scale_outputs <- foreach(unscaled_layer=unstack(x), 
                                 .packages=c('raster', 'teamlucc'),
                                 .export=c('scale_layer')) %dopar% {
            scale_output <- scale_layer(unscaled_layer, power_of, max_out, 
                                        round_output, do_scaling)
        }
    }
    names(scale_outputs) <- names(x)
    return(scale_outputs)
}

#' @rdname scale_raster
#' @aliases scale_raster,RasterStack,ANY-method
setMethod("scale_raster", signature(x="RasterStack"),
    function(x, power_of, max_out, round_output, do_scaling) {
        ret <- scale_stack_or_brick(x, power_of, max_out, round_output, 
                                    do_scaling)
        return(ret)
    }
)

#' @rdname scale_raster
#' @aliases scale_raster,RasterBrick,ANY-method
setMethod("scale_raster", signature(x="RasterBrick"),
    function(x, power_of, max_out, round_output, do_scaling) {
        ret <- scale_stack_or_brick(x, power_of, max_out, round_output, 
                                    do_scaling)
        return(ret)
    }
)
