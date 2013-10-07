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
#' its relationships with biomass estimation in Rondônia, Brazilian Amazon.  
#' Acta Amazonica 35:249-257.
#'
#' Gonzalez, R. C. 2008. Digital image processing. 3rd ed. Prentice Hall, Upper 
#' Saddle River, N.J, pages 830-836.
#'
#' Haralick, R. M., K. Shanmugam, and I. Dinstein. 1973. Textural features for 
#' image classification. IEEE Transactions on Systems, Man and Cybernetics 
#' SMC-3:610-621.
#'
#' Pratt, W. K. 2007. Digital image processing: PIKS Scientific inside4th ed., 
#' Newly updated and rev. ed. Wiley-Interscience, Hoboken, N.J pages 540-541, 
#' 563-566.
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
    x_grey <- cut(x, breaks=seq(cellStats(x, 'min'), cellStats(x, 'max'), 
                                length.out=n_grey + 1), include.lowest=TRUE)

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

    calc_texture <- function(rast, statistics, base_indices, offset_indices, n_grey, 
                             ...) {
        G <- matrix(0, n_grey, n_grey)
        co_occur <- cbind(rast[,,1][base_indices], rast[,,1][offset_indices])
        for (n in 1:nrow(co_occur)) {
            G[co_occur[n, 1], co_occur[n, 2]] <- G[co_occur[n, 1],
                                                   co_occur[n, 2]] + 1
        }
        Pij = G / sum(G)

        # Calculate rowSums and colSums of Pij
        rowsum = rowSums(Pij)
        colsum = colSums(Pij)
        # Calcuate mr and mc (forms of col and row means) and sig2r and sig2c 
        # (measures of row and column variance)
        mr = sum(c(1:nrow(Pij)) * rowsum)
        mc = sum(c(1:ncol(Pij)) * colsum)
        sig2r = sum((c(1:nrow(Pij)) - mr)^2 * rowsum)
        sig2c = sum((c(1:ncol(Pij)) - mc)^2 * colsum)

        # Make a matrix of i's and a matrix of j's to be used in the below 
        # matrix calculations. These matrices are the same shape as Pij with 
        # the entries equal to the i indices of each cell (for the imat matrix, 
        # which is indexed over the rows) or the j indices of each cell (for 
        # the jmat matrix, which is indexed over the columns).
        imat <- matrix(rep(1:nrow(Pij), ncol(Pij)), nrow=nrow(Pij))
        jmat <- matrix(rep(1:ncol(Pij), nrow(Pij)), ncol=ncol(Pij), byrow=TRUE)

        textures <- c()
        if ('mean' %in% statistics) {
            # Defined as in Lu and Batistella, 2005, page 252
            textures <- c(textures, mr)
        }
        if ('variance' %in% statistics) {
            # Defined as in Haralick, 1973, page 619 (equation 4)
            textures <- c(textures, sum(rowSums((imat - mr)^2 * Pij)))
        }
        if ('covariance' %in% statistics) {
            # Defined as in Pratt, 2007, page 540
            textures <- c(textures, sum(rowSums((imat - mr) *
                                                (jmat - mc) * Pij)))
        }
        if ('homogeneity' %in% statistics) {
            # Defined as in Gonzalez and Woods, 2009, page 832
            textures <- c(textures, sum(rowSums(Pij / (1 + abs(imat - jmat)))))
        }
        if ('contrast' %in% statistics) {
            # Defined as in Gonzalez and Woods, 2009, page 832
            textures <- c(textures, sum(rowSums((imat - jmat)^2 * Pij)))
        }
        if ('dissimilarity' %in% statistics) {
            #TODO: Find source for dissimilarity
            stop('Error: dissimilarity not yet supported')
        }
        if ('entropy' %in% statistics) {
            # Defined as in Gonzalez and Woods, 2009, page 832
            textures <- c(textures, -sum(rowSums(Pij * log2(Pij + .001))))
        }
        if ('second_moment' %in% statistics) {
            # Defined as in Haralick, 1973, page 619
            textures <- c(textures, sum(rowSums(Pij^2)))
        }
        if ('correlation' %in% statistics) {
            # Defined as in Gonzalez and Woods, 2009, page 832
            textures <- c(textures, sum(rowSums(((imat - mr) * (jmat - mc) * Pij) /
                                                (sqrt(sig2r) * sqrt(sig2c)))))
        }
    }

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
