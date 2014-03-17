#' R/C++ implementation of Xiaolin Zhu's CLOUD_REMOVE IDL script
#'
#' @export
#' @import Rcpp
#' @import RcppArmadillo
#' @importFrom spatial.tools rasterEngine sfQuickInit sfQuickStop
#' @param cloudy the cloudy image (base image) as a \code{Raster*}
#' @param clear the clear image as a \code{Raster*} to use for filling 
#' \code{img_cloudy}
#' @param cloud_mask the cloud mask as a \code{RasterLayer}, with each cloud 
#' patch assigned a unique integer code. Areas that are clear in both 
#' \code{cloudy_rast} and \code{clear_rast} should be coded 0, while areas that 
#' are clouded in \code{clear_rast} should be coded -1.
#' @param out_name filename for cloud filled image
#' @param num_class set the estimated number of classes in image
#' @param min_pixel the sample size of similar pixels
#' @param max_pixel the maximum sample size to search for similar pixels
#' @param cloud_nbh the range of cloud neighborhood (in pixels)
#' @param DN_min the minimum valid DN value
#' @param DN_max the maximum valid DN value
#' @param ... additional arguments to pass to rasterEngine
#' @references Zhu, X., Gao, F., Liu, D., Chen, J., 2012. A modified
#' neighborhood similar pixel interpolator approach for removing thick clouds 
#' in Landsat images. Geoscience and Remote Sensing Letters, IEEE 9, 521--525.
cloud_remove <- function(cloudy, clear, cloud_mask, out_name, num_class=1, 
                         min_pixel=20, max_pixel=1000, cloud_nbh=1, DN_min=0, 
                         DN_max=255, ...) {
    if (nlayers(cloudy) != nlayers(clear)) {
        stop('number of layers in cloudy_rast must match number of layers in clear_rast')
    }
    if (nlayers(cloud_mask) != 1) {
        stop('mask_rast should have only one layer')
    }
    
    # Wrapper around C++ cloud fill function
    cloud_fill_wrapper <- function(cloudy, clear, cloud_mask, num_class, 
                                   min_pixel, max_pixel, cloud_nbh, DN_min, 
                                   DN_max, ...) {
        dims=dim(cloudy)
        # RcppArmadillo crashes when you pass it a cube, so resize and pass 
        # mats
        cloudy <- array(cloudy, dim=c(dims[1] * dims[2], dims[3]))
        clear <- array(clear, dim=c(dims[1] * dims[2], dims[3]))
        cloud_mask <- array(cloud_mask, dim=c(dims[1] * dims[2]))
        filled <- cloud_fill(cloudy, clear, cloud_mask, dims, num_class, 
                             min_pixel, max_pixel, cloud_nbh, DN_min, DN_max)
        # RcppArmadillo crashes when you return a cube, so resize the returned 
        # mat
        filled <- array(filled, dim=c(dims[1], dims[2], dims[3]))
    }

    # bs <- blockSize(cloudy)
    # out <- brick(cloudy, values=FALSE)
    # out <- writeStart(out, rasterTmpFile())
    # message(paste0("**", bs$n, " blocks to process**"))
    # for (block_num in 1:bs$n) {
    #     message("Processing block ", block_num, " - ", appendLF=FALSE)
    #     cloudy_bl <- getValuesBlock(cloudy, row=bs$row[block_num], nrows=bs$nrows[block_num])
    #     clear_bl <- getValuesBlock(clear, row=bs$row[block_num], nrows=bs$nrows[block_num])
    #     cloud_mask_bl <- getValuesBlock(cloud_mask, row=bs$row[block_num], nrows=bs$nrows[block_num])

    #     dims <- c(bs$nrows[block_num], ncol(cloudy), nlayers(cloudy))

    #     filled <- cloud_fill(cloudy_bl, clear_bl, cloud_mask_bl, dims, 
    #                          num_class, min_pixel, max_pixel, cloud_nbh, 
    #                          DN_min, DN_max)

    #     out <- writeValues(out, filled, bs$row[block_num])
    # }
    # out <- writeStop(out)

    out <- rasterEngine(cloudy=cloudy, clear=clear, 
                        cloud_mask=cloud_mask,
                        fun=cloud_fill_wrapper,
                        args=list(num_class=num_class, min_pixel=min_pixel, 
                        max_pixel=max_pixel, cloud_nbh=cloud_nbh, 
                        DN_min=DN_min, DN_max=DN_max),
                        processing_unit='chunk',
                        outbands=nlayers(cloudy), outfiles=1, ...)

    return(out)
}
