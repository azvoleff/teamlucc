#' Image texture measures from grey-level co-occurrence matrices (GLCM)
#'
#' @export
#' @param x a /code{RasterLayer}
#' @param n_grey number of grey levels to use in texture calculation
#' @param window the window size to consider for texture calculation as a two 
#' element integer vector
#' @param shift a 2 element integer vector giving the shift (Q in Gonzalez and 
#' Woods, 2008). 
#' @param statistics A list of GLCM texture measures to calculate.  Can include 
#' any (one or more) of the following: 'mean', 'variance', 'covariance', 
#' 'homogeneity', 'contrast', 'dissimilarity', 'entropy', 'second_moment', 
#' 'correlation'.
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
#' \dontrun{
#' L5TSR_1986 <- stack(system.file('extdata/L5TSR_1986.dat', package='teamr'))
#'
#' x <- raster(L5TSR_1986, layer=1)
#' textures <- glcm(x)
#' }
glcm <- function(x, n_grey=32, window=c(3, 3), shift=c(1, 1),
                 statistics=c('mean', 'variance', 'covariance', 'homogeneity', 
                              'contrast', 'dissimilarity', 'entropy', 
                              'second_moment', 'correlation')) {
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

    # Resample the image to the required number of grey levels
    message(paste('Resampling to', n_grey, 'grey levels'))
    x_grey <- cut(x, breaks=seq(cellStats(x, 'min'), cellStats(x, 'max'), 
                                length.out=n_grey + 1), include.lowest=TRUE)

    message('Calculating textures')
    # Calculate column major indices for base and offset images
    ul_base <- c(1, 1)
    if (shift[1] < 0) {ul_base[1] <- ul_base[1] + abs(shift[1])}
    if (shift[2] < 0) {ul_base[2] <- ul_base[2] + abs(shift[2])}
    ul_offset <- ul_base + shift
    base_cols <- matrix(rep(ul_base[2] + seq(1, window[2]), window[1]) - 1, 
                        nrow=window[1], byrow=TRUE)
    base_rows <- matrix(rep(ul_base[1] + seq(1, window[1]), window[2]) - 1, 
                        nrow=window[1])
    offset_cols <- matrix(rep(ul_offset[2] + seq(1, window[2]), window[1]) - 1, 
                          nrow=window[1], byrow=TRUE)
    offset_rows <- matrix(rep(ul_offset[1] + seq(1, window[1]), window[2]) - 1, 
                          nrow=window[1])
    offset_indices <- (offset_cols - 1) * (window[1]+abs(shift[1])) + offset_rows
    base_indices <- (base_cols - 1) * (window[1]+abs(shift[1])) + base_rows

    texture_img <- rasterEngine(rast=x_grey, fun=calc_texture,
                                args=list(statistics=statistics, 
                                          base_indices=base_indices, 
                                          offset_indices=offset_indices,
                                          n_grey=n_grey), 
                                window_dims=(window + abs(shift)))

    names(texture_img) <- statistics 
    texture_img <- setMinMax(texture_img)

    return(texture_img)
}
n_grey=32
window=c(3, 3)
shift=c(1, 1)
statistics=c('mean', 'variance', 'covariance', 'homogeneity', 
             'contrast', 'dissimilarity', 'entropy', 
             'second_moment', 'correlation')
L5TSR_1986 <- stack(system.file('extdata/L5TSR_1986.dat', package='teamr'))
x <- raster(L5TSR_1986, layer=1)

textures <- glcm(x)
