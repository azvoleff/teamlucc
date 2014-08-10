#' Change Direction Image for CVAPS
#'
#' This code calculate the change direction image for the Change Vector 
#' Analysis in Posterior Probability Space (CVAPS) method of Chen et al. 2011.  
#' Use the change direction image in conjunction with the change magnitude 
#' image from \code{chg_dir}, and \code{DFPS} to use the Double Window Flexible 
#' Pace Search method (Chen et al. 2003) to determine the threshold to use to 
#' map areas of change and no-change.
#'
#' @export
#' @importFrom spatial.tools rasterEngine
#' @importFrom tools file_path_sans_ext
#' @param t1p time 0 posterior probability \code{Raster*}
#' @param t2p time 1 posterior probability \code{Raster*}
#' @param filename (optional) filename for output change direction
#' \code{RasterLayer}
#' @param overwrite whether to overwrite existing files (otherwise an error 
#' will be raised)
#' @param verbose whether to print detailed status messages
#' @param ... additional parameters to pass to rasterEngine
#' @return \code{Raster*} object with change direction image
#' @references Chen, J., P. Gong, C. He, R. Pu, and P. Shi. 2003.
#' Land-use/land-cover change detection using improved change-vector analysis.
#' Photogrammetric Engineering and Remote Sensing 69:369-380.
#' 
#' Chen, J., X. Chen, X. Cui, and J. Chen. 2011. Change vector analysis in 
#' posterior probability space: a new method for land cover change detection.  
#' IEEE Geoscience and Remote Sensing Letters 8:317-321.
#' @examples
#' \dontrun{
#' t0_train_data <- get_pixels(L5TSR_1986, L5TSR_1986_2001_training, "class_1986",training=.6)
#' t0_model <- train_classifier(t0_train_data)
#' t0_preds <- classify(L5TSR_1986, t0_model)
#' t1_train_data <- get_pixels(L5TSR_2001, L5TSR_1986_2001_training, "class_2001", training=.6)
#' t1_model <- train_classifier(t1_train_data)
#' t1_preds <- classify(L5TSR_2001, t1_model)
#' t0_t1_chgdir <- chg_dir(t0_preds$probs, t1_preds$probs)
#' }
chg_dir <- function(t1p, t2p, filename, overwrite=FALSE, verbose=FALSE, ...) {
    if (proj4string(t1p) != proj4string(t2p)) {
        stop('t0 and t1 coordinate systems do not match')
    }
    if (extent(t1p) != extent(t2p)) {
        stop('t0 and t1 extents do not match')
    }
    if (nlayers(t1p) != nlayers(t2p)) {
        stop('t0 and t1 probability maps have differing number of classes')
    }
    if (!missing(filename) && file_test('-f', filename) && !overwrite) {
        stop('output file already exists and overwrite=FALSE')
    }

    n_classes <- nlayers(t1p)
    if (n_classes == 1) {
        stop('cannot calculate change probabilities for only one class')
    }

    # calc_chg_dir_st <- function(t1p, t2p, n_classes, ...) {
    #     # Calculate change direction (eqns 5 and 6 in Chen 2011)
    #     dP <- array(t2p - t1p, dim=c(dim(t1p)[1], dim(t1p)[2], n_classes))
    #     unit_vecs <- array(diag(n_classes), dim=c(n_classes, n_classes))
    #     Eab <- apply(dP, c(1, 2), function(pixel) pixel %*% unit_vecs)
    #     chgdir <- apply(Eab, c(2, 3),
    #                     function(pixel) which(pixel == max(pixel)))
    #     chgdir <- array(chgdir, dim=c(dim(t1p)[1], dim(t1p)[2], 1))
    #     return(chgdir)
    # }
    # out <- rasterEngine(t1p=t1p, t2p=t2p, fun=calc_chg_dir_st, 
    #                     args=list(n_classes=n_classes), 
    #                     outbands=1, datatype='INT2S', ...)
    #
    # # spatial.tools can only output the raster package grid format - so output 
    # # to a tempfile in that format then copy over to the requested final output 
    # # format if a filename was supplied
    # if (!missing(filename)) {
    #     out <- writeRaster(out, filename=filename, dataType='INT2S', 
    #                        overwrite=overwrite)
    # }
    
    if (missing(filename)) {
        filename <- rasterTmpFile()
        overwrite <- TRUE
    }
   
    bs <- blockSize(t1p)
    out <- raster(t1p)
    out <- writeStart(out, filename=filename, overwrite=overwrite)
    for (block_num in 1:bs$n) {
        if (verbose > 0) {
            message("Processing block ", block_num, " of ", bs$n, "...")
        }
        dims <- c(bs$nrows[block_num], ncol(t1p), nlayers(t1p))
        t1p_bl <- array(getValuesBlock(t1p, row=bs$row[block_num], 
                                 nrows=bs$nrows[block_num]),
                        dim=c(dims[1] * dims[2], dims[3]))
        t2p_bl <- array(getValuesBlock(t2p, row=bs$row[block_num], 
                                       nrows=bs$nrows[block_num]),
                        dim=c(dims[1] * dims[2], dims[3]))
        chg_dirs <- calc_chg_dir(t1p_bl, t2p_bl)
        out <- writeValues(out, chg_dirs, bs$row[block_num])
    }
    out <- writeStop(out)

    return(out)
}
