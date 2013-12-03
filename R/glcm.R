#' Image texture measures from grey-level co-occurrence matrices (GLCM)
#'
#' @export
#' @import spatial.tools Rcpp RcppArmadillo
#' @param layer a /code{RasterLayer}
#' @param n_grey number of grey levels to use in texture calculation
#' @param window the window size to consider for texture calculation as a two 
#' element integer vector
#' @param shift a 2 element integer vector giving the shift (Q in Gonzalez and 
#' Woods, 2008). 
#' @param statistics A list of GLCM texture measures to calculate.  Can include 
#' any (one or more) of the following: 'mean', 'mean_ENVI', 'variance', 
#' 'variance_ENVI', 'homogeneity', 'contrast', 'dissimilarity', 'entropy', 
#' 'second_moment', and/or 'correlation'.
#' @param ... additional parameters to pass to rasterEngine
#' @return A /code{RasterLayer} with the calculated requested GLCM texture 
#' measures.
#' @references Lu, D., and M. Batistella. 2005. Exploring TM image texture and 
#' its relationships with biomass estimation in Rond\^{o}nia, Brazilian Amazon.  
#' Acta Amazonica 35:249-257.
#'
#' Gonzalez, R. C. 2008. Digital image processing. 3rd ed. Prentice Hall, Upper 
#' Saddle River, N.J, pages 830-836.
#'
#' Haralick, R. M., K. Shanmugam, and I. Dinstein. 1973. Textural features for 
#' image classification. IEEE Transactions on Systems, Man and Cybernetics 
#' SMC-3:610-621.
#'
#' Pratt, W. K. 2007. Digital image processing: PIKS Scientific inside. 4th ed.
#' Wiley-Interscience, Hoboken, N.J pages 540-541, 563-566.
#' @examples
#' L5TSR_1986 <- stack(system.file('extdata/L5TSR_1986.dat', package='teamr'))
#' textures <- glcm(raster(L5TSR_1986, layer=1))
#' plot(textures)
glcm <- function(layer, n_grey=32, window=c(3, 3), shift=c(1, 1),
                 statistics=c('mean', 'variance', 'homogeneity', 'contrast', 
                              'dissimilarity', 'entropy', 'second_moment', 
                              'correlation'), ...) {
    if (length(window) != 2) {
        stop('window must be integer vector of length 2')
    }
    if (length(shift) != 2) {
        stop('shift must be integer vector of length 2')
    }
    if ((window[1] < 3) || (window[2] < 3)) {
        stop('both elements of window must be  >= 3')
    }
    if ((window[1] %% 2 == 0) || (window[2] %% 2 == 0)) {
        stop('both elements of window must be odd')
    }
    avail_stats <- c('mean', 'mean_ENVI', 'variance', 'variance_ENVI', 
                     'homogeneity', 'contrast', 'dissimilarity', 'entropy', 
                     'second_moment', 'correlation')
    stat_check <- unlist(lapply(statistics, function(x) x %in% avail_stats))
    if (sum(stat_check) != length(stat_check)) {
        stop(paste('invalid texture(s):',
                   paste(statistics[!stat_check], collapse=', ')))
    }

    # Resample the image to the required number of grey levels
    #message(paste('Resampling to', n_grey, 'grey levels...'))
    layer <- raster::cut(layer, breaks=seq(cellStats(layer, 'min'), 
                                           cellStats(layer, 'max'), 
                                           length.out=n_grey + 1), 
                         include.lowest=TRUE)

    #message('Calculating textures...')
    texture_img <- calc_texture_full_image(raster::as.matrix(layer), 
                                           n_grey, window, shift, statistics)
    texture_img <- stack(apply(texture_img, 3, raster, template=layer))

    names(texture_img) <- statistics 

    return(texture_img)
}
