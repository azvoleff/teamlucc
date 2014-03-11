#' Xiaolin Zhu's CLOUD_REMOVE.pro script implemented in R
#'
#' @importFrom spatial.tools rasterEngine
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
#' @param cloud_nbh the range of cloud neighborhood (in pixels)
#' @param DN_min the minimum valid DN value
#' @param DN_max the maximum valid DN value
#' @references Zhu, X., Gao, F., Liu, D., Chen, J., 2012. A modified
#' neighborhood similar pixel interpolator approach for removing thick clouds 
#' in Landsat images. Geoscience and Remote Sensing Letters, IEEE 9, 521--525.
cloud_remove <- function(cloudy, clear, cloud_mask, out_name, 
                         num_class=1, min_pixel=20, cloud_nbh=1, DN_min=0, 
                         DN_max=255) {

    if (nlayers(cloudy) != nlayers(clear)) {
        stop('number of layers in cloudy_rast must match number of layers in clear_rast')
    }
    if (nlayers(cloud_mask) != 1) {
        stop('mask_rast should have only one layer')
    }

    # # Loop over blocks
    # bs <- blockSize(cloudy_rast)
    # block_num <- 1
    # cloudy_bl <- getValues(cloudy_rast, row=bs$row[block_num], nrows=bs$nrows[block_num])
    # clear_bl <- getValues(clear_rast, row=bs$row[block_num], nrows=bs$nrows[block_num])
    # mask_bl <- getValues(mask_rast, row=bs$row[block_num], nrows=bs$nrows[block_num])
    # Note that the first dimension of the array is columns (of original 
    # image), second is rows (of original image), and third is layers
    # cloudy_bl <- array(cloudy_bl, c(ncol(cloudy_rast), bs$nrows[block_num], nlayers(cloudy_rast)))
    # clear_bl <- array(clear_bl, c(ncol(clear_rast), bs$nrows[block_num], nlayers(clear_rast)))
    # mask_bl <- array(mask_bl, c(ncol(mask_rast), bs$nrows[block_num]))
    # image(cloudy_bl[,,1])
    # image(clear_bl[,,1])
    # image(mask_bl)

    cloud_remove_block <- function(cloudy, clear, cloud_mask, num_class=1, min_pixel=20, 
                           cloud_nbh=1, DN_min=0, DN_max=255, ...) {
        num_bands <- dim(cloudy)[3]
        # Loop over clouds
        cloud_codes <- unique(as.vector(cloud_mask))
        cloud_codes <- cloud_codes[cloud_codes != 0] # 0 is code for background
        for (cloud_code in cloud_codes) {
            cloud_indices <- which(cloud_mask == cloud_code, arr.ind=TRUE)

            # Calculate the neighborhood of cloud Here the coordinates are all 
            # relative to the block, NOT the raster. Remember the blocks are 
            # flipped over the horizontal plane compared to the original 
            # raster. In other words, row 1, column  1 of the raster (upper 
            # left) is located at 1, 1 in Cartesian coordinates (lower left) 
            # when plotted, and is at row nrow(block), column 1, of the matrix.
            #
            # Also remember that image() plots the image rotated 
            # counterclockwise 90 degrees relative to the normal printed 
            # version of a matrix.
            left_cloud <- min(cloud_indices[, 2])
            right_cloud <- max(cloud_indices[, 2])
            up_cloud <- min(cloud_indices[, 1])
            down_cloud <- max(cloud_indices[, 1])

            # Now add in the cloud neighborhood
            left_region <- max(c(1, left_cloud - cloud_nbh))
            right_region <- min(c(ncol(cloudy), right_cloud + cloud_nbh))
            up_region <- max(c(1, up_cloud - cloud_nbh))
            down_region <- min(c(nrow(cloudy), down_cloud + cloud_nbh))

            sub_cloudy <- cloudy[up_region:down_region, left_region:right_region, ]
            sub_clear <- clear[up_region:down_region, left_region:right_region, ]
            sub_cloud_mask <- cloud_mask[up_region:down_region, left_region:right_region, ]

            # Calculate position of cloud center
            x_center <- ncol(sub_cloudy)/ 2.0
            y_center <- nrow(sub_cloudy)/ 2.0

            # image(sub_cloudy[, , 1])
            # image(sub_clear[, , 1])
            # image(sub_cloud_mask)

            # Compute the threshold for what is a "similar" pixel
            similar_th_band <- apply(sub_clear, 3, sd)
            similar_th_band <- similar_th_band * 2 / num_class

            # Find indices of clear and clouded pixels
            is_cloud <- sub_cloud_mask == cloud_code
            is_clear <- sub_cloud_mask == 0

            sub_cloud_mask_clear_indices <- which(is_clear, arr.ind=TRUE)
            num_clear_pixels <- nrow(sub_cloud_mask_clear_indices)
            sub_cloud_mask_cloud_indices <- which(is_cloud, arr.ind=TRUE)
            num_cloud_pixels <- nrow(sub_cloud_mask_cloud_indices)

            sub_clear_clear <- apply(sub_clear, 3, function(x) x[is_clear])
            sub_cloudy_clear <- apply(sub_cloudy, 3, function(x) x[is_clear])

            for (ic in 1:num_cloud_pixels) {
                # Calculate row and column location of target pixel
                ri <- sub_cloud_mask_cloud_indices[ic, 1]
                ci <- sub_cloud_mask_cloud_indices[ic, 2]

                # Calculate distance between target pixel and center of cloud
                r2 <- ((x_center - ri)^2 + (y_center - ci)^2)^0.5
                clear_dists <- ((sub_cloud_mask_clear_indices[, 1] - ri)^2 +
                                (sub_cloud_mask_clear_indices[, 2] - ci)^2)^0.5

                order_clear <- order(clear_dists)
                # Below line avoids comparing a pixel with itself
                order_clear <- order_clear[-1]

                # Below needs to be coded in C++ for speed
                iclear <- 1
                num_similar <- 0
                cloudy_similar <- matrix(0, min_pixel, num_bands)
                clear_similar <- matrix(0, min_pixel, num_bands)
                rmse_similar <- matrix(0, min_pixel) # Based on spectral distance
                dis_similar <- matrix(0, min_pixel) # Based on spatial distance
                while(num_similar <= (min_pixel-1) && iclear <= (num_clear_pixels-1)) {
                    # In below line: sub_clear[ri, ci, iband] is the clear 
                    # image value of the clouded pixel we are focusing on
                    indicate_similar <- sum((sub_clear_clear[order_clear[iclear], ]
                                             - sub_clear[ri, ci, ]) <= similar_th_band)
                    # Below only runs if there are similar pixels in all bands
                    if (indicate_similar == num_bands) {
                        cloudy_similar[num_similar + 1, ] <- sub_cloudy_clear[order_clear[iclear], ]
                        clear_similar[num_similar + 1, ] <- sub_clear_clear[order_clear[iclear], ]
                        rmse_similar[num_similar + 1] <- (sum((sub_clear_clear[order_clear[iclear], ] - sub_clear[ri, ci, ])^2) / num_bands)^0.5
                        dis_similar[num_similar + 1] <- clear_dists[order_clear[iclear]]
                        num_similar <- num_similar + 1
                    }
                    iclear <- iclear + 1
                }

                if (num_similar > 1) {
                    rmse_similar <- rmse_similar[1:num_similar, ]
                    dis_similar <- dis_similar[1:num_similar, ]
                    cloudy_similar <- cloudy_similar[1:num_similar, ]
                    clear_similar <- clear_similar[1:num_similar, ]
                    rmse_similar_norm <- (rmse_similar - min(rmse_similar)) / (max(rmse_similar) - min(rmse_similar) + 0.000001) + 1.0
                    dis_similar_norm <- (dis_similar - min(dis_similar)) / (max(dis_similar) - min(dis_similar) + 0.000001) + 1.0
                    C_D <- (rmse_similar_norm) * dis_similar_norm + 0.0000001
                    weight <- (1.0 / C_D) / sum(1.0 / C_D)

                    # Compute the time weight
                    W_T1 <- r2 / (r2 + mean(dis_similar))
                    W_T2 <- mean(dis_similar) / (r2 + mean(dis_similar))

                    # Make predictions
                    predict_1 <- apply(cloudy_similar * weight, 2, sum)
                    predict_2 <- sub_clear[ri, ci, ] + apply((cloudy_similar - clear_similar)*weight, 2, sum)
                    for (iband in 1:num_bands) {
                        if (predict_2[iband] > DN_min && predict_2[iband] < DN_max) {
                            sub_cloudy[ri, ci, iband] <- W_T1 * predict_1[iband] + W_T2 * predict_2[iband]
                        } else {
                            sub_cloudy[ri, ci, iband] <- predict_1[iband]
                        }
                    }

                } else {
                    # If no similar pixel, use mean of all pixels in cloud 
                    # neighborhood for a simple linear adjustment
                    mean_diff <- apply(sub_cloudy_clear - sub_clear_clear, 2, mean)
                    sub_cloudy[ri, ci, ] <- sub_clear[ri, ci, ] + mean_diff
                }
            }
            cloudy[up_region:down_region, left_region:right_region, ] <- sub_cloudy
        }
        return(cloudy)
    }

    out <- rasterEngine(cloudy=cloudy, clear=clear, 
                        cloud_mask=cloud_mask,
                        fun=cloud_remove_block,
                        processing_unit='chunk',
                        args=list(num_class=num_class, min_pixel=min_pixel, 
                                  cloud_nbh=cloud_nbh, DN_min=DN_min, 
                                  DN_max=DN_max),
                        filename=out_name,
                        outbands=nlayers(cloudy), verbose=TRUE)

    return(out)
}
