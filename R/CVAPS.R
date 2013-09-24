#' Change Vector Analysis in Posterior Probability Space
#'
#' This code implements the Change Vector Analysis in Posterior Probability 
#' Space method of Chen et al. 2011. This function produces the change 
#' magnitude and change direction images. Use the change magnitude image in 
#' conjunction with the time 0 and time 1 images to produce a shapefile of 
#' change polygons, then run \code{DFPS} to use the Double Window Flexible Pace 
#' Search method (Chen et al. 2003) to determine the threshold to use to map 
#' areas of change and no-change.
#'
#' @export
#' @param t0_base base name for time 0 classification and posterior probability 
#' files
#' @param t1_base base name for time 1 classification and posterior probability 
#' files
#' @param out_file_base the base name to use when naming the output files
#' @references Chen, J., P. Gong, C. He, R. Pu, and P. Shi. 2003.
#' Land-use/land-cover change detection using improved change-vector analysis.
#' Photogrammetric Engineering and Remote Sensing 69:369-380.
#' 
#' Chen, J., X. Chen, X. Cui, and J. Chen. 2011. Change vector analysis in 
#' posterior probability space: a new method for land cover change detection.  
#' IEEE Geoscience and Remote Sensing Letters 8:317-321.
CVAPS <- function(t0_base, t1_base, out_file_base) {
    t0_prob <- brick(paste(t0_base, '_probs.envi', sep=''))
    t1_prob <- brick(paste(t1_base, '_probs.envi', sep=''))

    if (proj4string(t0_prob) != proj4string(t1_prob)) {
        stop('Error: t0 and t1 coordinate systems do not match')
    }
    if (extent(t0_prob) != extent(t1_prob)) {
        stop('Error: t0 and t1 extents do not match')
    }
    if (nlayers(t0_prob) != nlayers(t1_prob)) {
        stop('Error: t0 and t1 probability maps have differing number of classes')
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
