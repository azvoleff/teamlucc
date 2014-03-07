#' @import raster
#' @import spatial.tools
#' @param cloudy_rast the cloudy image (base image) as a \code{Raster*}
#' @param clear_rast the clear image as a \code{Raster*} to use for filling 
#' \code{img_cloudy}
#' @param mask_rast the cloud mask as a \code{RasterLayer}, with each cloud 
#' patch assigned a unique integer code. Areas that are clear in both 
#' \code{cloudy_rast} and \code{clear_rast} should be coded 0, while areas that 
#' are clouded in \code{clear_rast} should be coded -1.
#' @param num_class set the estimated number of classes in image
#' @param min_pixel the sample size of similar pixels
#' @param cloud_nbh the range of cloud neighborhood (in pixels)
#' @param DN_min minimum of DN
#' @param DN_max maximum of DN
cloud_remove <- function(cloudy_rast, clear_rast, mask_rast, out_name, 
                         num_class=1, min_pixel=20, cloud_nbh=1, DN_min=0, 
                         DN_max=255) {

    if (nlayers(cloudy_rast) != nlayers(clear_rast)) {
        stop('number of layers in cloudy_rast must match number of layers in clear_rast')
    }
    if (nlayers(mask_rast) != 1) {
        stop('mask_rast should have only one layer')
    }

    out <- brick(cloudy_rast)
    out <- writeStart(out, out_name)

    num_bands <- nlayers(cloudy_rast)

    # Loop over blocks
    bs <- blockSize(cloudy_rast)
    for (block_num in 1:bs$n) {
        cloudy_bl <- getValues(cloudy_rast, row=bs$row[block_num], nrows=bs$nrows[block_num])
        clear_bl <- getValues(clear_rast, row=bs$row[block_num], nrows=bs$nrows[block_num])
        mask_bl <- getValues(mask_rast, row=bs$row[block_num], nrows=bs$nrows[block_num])

        # Note that the first dimension of the array is columns (of original 
        # image), second is rows (of original image), and third is layers
        cloudy_bl <- array(cloudy_bl, c(ncol(cloudy_rast), bs$nrows[block_num], nlayers(cloudy_rast)))
        clear_bl <- array(clear_bl, c(ncol(clear_rast), bs$nrows[block_num], nlayers(clear_rast)))
        mask_bl <- array(mask_bl, c(ncol(mask_rast), bs$nrows[block_num], nlayers(mask_rast)))

        image(cloudy_bl[,,1])
        image(clear_bl[,,1])
        image(mask_bl[,,1])

        # Loop over clouds (make more efficient by replacing below with a 
        # call to unique and by only looping over codes that actually occur 
        # in this patch
        cloud_codes <- unique(as.vector(mask_bl))
        cloud_codes <- cloud_codes[cloud_codes != 0] # 0 is code for background
        for (cloud_code in cloud_codes) {
            cloud_indices <- which(mask_bl == cloud_code, arr.ind=TRUE)

            left_cloud <- min(cloud_indices[, 1])
            right_cloud <- max(cloud_indices[, 1])
            up_cloud <- min(cloud_indices[, 2]) # remember row1 is on top
            down_cloud <- max(cloud_indices[, 2])

            # Neighborhood of cloud - again remember nrow(cloudy_bl) is equal 
            # to ncol(cloudy_rast)
            left_region <- max(c(1, left_cloud - cloud_nbh + 1))
            right_region <- min(c(nrow(cloudy_bl), right_cloud + cloud_nbh))
            up_region <- max(c(1, up_cloud - cloud_nbh + 1))
            down_region <- min(c(ncol(cloudy_bl), down_cloud + cloud_nbh))

            a_region <- right_region - left_region + 1
            b_region <- down_region - up_region + 1
            # Calculate position of cloud center
            x_center <- a_region / 2.0
            y_center <- b_region / 2.0

            sub_cloudy_bl <- cloudy_bl[left_region:right_region, up_region:down_region, ]
            sub_clear_bl <- clear_bl[left_region:right_region, up_region:down_region, ]
            sub_mask_bl <- mask_bl[left_region:right_region, up_region:down_region, ]

            # Compute the threshold for what is a "similar" pixel
            similar_th_band <- apply(sub_clear_bl, 3, sd)
            similar_th_band <- similar_th_band * 2 / num_class

            # Find indices of clear and clouded pixels
            is_cloud <- sub_mask_bl == cloud_code
            is_clear <- sub_mask_bl == 0

            clear_indices <- which(is_clear, arr.ind=TRUE)
            num_clear_pixels <- nrow(clear_indices)
            cloud_indices <- which(is_cloud, arr.ind=TRUE)
            num_cloud_pixels <- nrow(cloud_indices)

            sub_clear_bl_clear <- sub_clear_bl[is_clear]
            sub_cloudy_bl_clear <- sub_cloudy_bl[is_clear]

            for (ic in 1:num_cloud_pixels) {
                # Calculate row and column location of target pixel
                ri <- cloud_indices[ic, 1]
                ci <- cloud_indices[ic, 2]

                # Calculate distance between target pixel and center of cloud
                r2 <- ((x_center - ci)^2 + (y_center - ri)^2)^0.5
                s_dis <- ((clear_indices[, 1] - ci)^2+(clear_indices[, 2] - ri)^2)^0.5

                order_clear <- order(s_dis)

                # Below needs to be coded in C++ for speed
                iclear <- 1
                isimilar <- 0
                cloudy_bl_similar <- matrix(0, min_pixel, num_bands)
                clear_bl_similar <- matrix(0, min_pixel, num_bands)
                rmse_similar <- matrix(0, min_pixel) # Based on spectral distance
                dis_similar <- matrix(0, min_pixel) # Based on spatial distance
                while (isimilar <= (min_pixel-1) && iclear <= (num_clear_pixels-1)) {
                    # Below line avoids comparing a pixel with itself
                    if (all(c(ci, ri) == clear_indices[order_clear[iclear], ])) next
                    # In below line: sub_clear_bl[ci, ri, iband] is the clear 
                    # image value of the clouded pixel we are focusing on
                    indicate_similar <- sum((sub_clear_bl_clear[order_clear[iclear]] - sub_clear_bl[ci, ri, ]) <= similar_th_band)
                    # Below only runs if there are similar pixels in all bands
                    #####################
                    #####################
                    #####################
                    #####################
                    # Got to here
                    #####################
                    #####################
                    if (indicate_similar == num_bands) {
                        for (iband in 1:num_bands) {
                            cloudy_bl_similar[isimilar,iband] <- sub_cloudy_bl_clear[, iband][order_clear[iclear]]
                            clear_bl_similar[isimilar,iband] <- sub_clear_bl_clear[, iband][order_clear[iclear]]
                        }
                        rmse_similar[isimilar] <- (sum((sub_clear_bl_clear[order_clear[iclear], ] - sub_clear_bl[ci,ri, ])^2) / float(num_bands))^0.5
                        dis_similar[isimilar] <- s_dis[order_clear[iclear]]
                        isimilar <- isimilar + 1
                    }
                    iclear <- iclear + 1
                }

                if (isimilar > 0) {
                    ind_null <- where(dis_similar > 0)
                    rmse_similar_norm <- (rmse_similar[ind_null] - min(rmse_similar[ind_null])) / (max(rmse_similar[ind_null]) - min(rmse_similar[ind_null]) + 0.000001) + 1.0
                    dis_similar_norm <- (dis_similar[ind_null] - min(dis_similar[ind_null])) / (max(dis_similar[ind_null]) - min(dis_similar[ind_null]) + 0.000001) + 1.0
                    C_D <- (rmse_similar_norm) * dis_similar_norm + 0.0000001
                    weight <- (1.0 / C_D) / sum(1.0 / C_D)

                    # Compute the time weight
                    W_T1 <- r2/(r2+mean(dis_similar[ind_null]))
                    W_T2 <- mean(dis_similar[ind_null])/(r2+mean(dis_similar[ind_null]))

                    # Make predictions
                    for (iband in 1:num_bands) {
                        predict_1 <- sum((cloudy_bl_similar[, iband])[ind_null] * weight)
                        predict_2 <- sub_clear_bl[ci, ri, iband] + sum(((cloudy_bl_similar[, iband])[ind_null]-(clear_bl_similar[, iband])[ind_null])*weight)
                        if (predict_2 > DN_min && predict_2 < DN_max) {
                            sub_cloudy_bl[ci, ri, iband] <- W_T1*predict_1+W_T2*predict_2
                        } else {
                            sub_cloudy_bl[ci, ri, iband] <- predict_1
                        }
                    }
                } else {
                    # If no similar pixel, use mean of all pixels in cloud 
                    # neighborhood for a simple linear adjustment
                    for (iband in 1:num_bands) {
                        mean_diff <- mean(sub_cloudy_bl_clear[, iband]-sub_clear_bl_clear[, iband])
                        sub_cloudy_bl[ci,ri,iband] <- sub_clear_bl[ci,ri,iband] + mean_diff
                    }
                }
            }
            # This needs to be done more carefully, to ensure only clouded 
            # pixels are replaced
            cloudy_bl[left_region:right_region,up_region:down_region, ] <- sub_cloudy_bl
        }
    }
    return(out)
}

