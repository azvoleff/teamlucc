#' Change Vector Analysis in Posterior Probability Space
#'
#' @export
#' @param t0_base base name for time 0 classification and posterior probability 
#' files
#' @param t1_base base name for time 1 classification and posterior probability 
#' files
#' @param out_file_base the base name to use when naming the output files
#' @references Chen, J., X. Chen, X. Cui, and J. Chen. 2011. Change vector 
#' analysis in posterior probability space: a new method for land cover change 
#' detection. IEEE Geoscience and Remote Sensing Letters 8:317-321.
CVAPS <- function(t0_base, t1_base, out_file_base) {
    t0_base <- 'C:/Users/azvoleff/SkyDrive/Publications/20140115 - TEAM LULC and Biodiversity/Code/L5TSR_1986'
    t1_base <- 'C:/Users/azvoleff/SkyDrive/Publications/20140115 - TEAM LULC and Biodiversity/Code/L5TSR_2001'
    out_file_base <- 'C:/Users/azvoleff/SkyDrive/Publications/20140115 - TEAM LULC and Biodiversity/Code/L5TSR_1986_to_2001'

    t0_prob <- brick(paste(t0_base, '_probs.envi', sep=''))
    t1_prob <- brick(paste(t1_base, '_probs.envi', sep=''))

    if (proj4string(t0_prob) != proj4string(t1_prob)) {
        stop('Error: t0 and t1 predicted class probability maps must have identical coordinate system')
    }
    if (nlayers(t0_prob) != nlayers(t1_prob)) {
        stop('Error: t0 and t1 predicted class probability maps must have identical number of classes')
    }
    
    n_classes <- nlayers(t0_prob)

    # Process over blocks to conserve memory. First setup the output rasters - 
    # one for the change magnitude image, the other for the change direction 
    # image.
    out_chgmag <- raster(t0_prob)
    out_chgmag <- writeStart(out_chgmag, paste(out_file_base, '_chgmag.envi', sep=''))
    out_chgdir <- raster(t0_prob)
    out_chgdir <- writeStart(out_chgdir, paste(out_file_base, '_chgdir.envi', sep=''))
    bs <- blockSize(t0_prob)
    for (block_num in 1:bs$n) {
        t0_prob_block <- getValues(t0_prob, row=bs$row[block_num], 
                                   nrows=bs$nrows[block_num])
        t1_prob_block <- getValues(t1_prob, row=bs$row[block_num], 
                                   nrows=bs$nrows[block_num])
        # Calculate change magnitude (eqn 3 in Chen 2011)
        dP <- t1_prob_block - t0_prob_block
        chgmag <- sqrt(rowSums((t1_prob_block - t0_prob_block)^2))
        writeValues(out_chgmag, chgmag, bs$row[block_num])
        # Calculate change vector (eqn 3 in Chen 2011). Use dot product of the 
        # base change vectors Ea,b,. Each base change vector represents 
        # complete change from class a to class b.
        Eab <- diag(n_classes)
        chgdir <- apply(dP %*% Eab, 1, function(r) which(r == max(r)))
        writeValues(out_chgdir, chgdir, bs$row[block_num])
    }
    out_chgmag <- writeStop(out_chgmag)
    out_chgdir <- writeStop(out_chgdir)
}
