context("cloud_remove")

# cloudy <- brick('C:/Users/azvoleff/Code/TEAM/teamlucc_testdata/xiaolinzhu/cloud remove update 20131023/test data/L20080724_cloudy')
# clear <- brick('C:/Users/azvoleff/Code/TEAM/teamlucc_testdata/xiaolinzhu/cloud remove update 20131023/test data/L20080606')
# cloud_mask <- raster('C:/Users/azvoleff/Code/TEAM/teamlucc_testdata/xiaolinzhu/cloud remove update 20131023/test data/cloud_mask')
# out_name <- 'C:/Users/azvoleff/Code/TEAM/teamlucc_testdata/xiaolinzhu/cloud remove update 20131023/test data/test_cloud_fill_in_R.envi'
#
# out_idl <- cloud_remove(cloudy, clear, cloud_mask, out_name)
#
# out_r <- cloud_remove(cloudy, clear, cloud_mask, out_name, use_IDL=FALSE)
#
# test_that("cloud_remove IDL and R results match", {
#           expect_equal(out_r, out_idl, tolerance=.001)
# })
#
# plot(cloudy[[1]])
# plot(clear[[1]])
# plot(out_idl[[1]])
# plot(out_r[[1]])
#
# image(cloudy_plt[, , 3])
# image(filled_plt[, , 3])
# image(clear_plt[, , 3])
# image(cloudy_plt[, , 3] - filled_plt[, , 3])
# image(filled_plt[, , 3] - clear_plt[, , 3])
# image(cloud_mask_plt)
