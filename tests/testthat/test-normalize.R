context("normalize")

suppressMessages(library(landsat))

###############################################################################
# Test on single band rasters
tl_res <- normalize(L5TSR_1986[[1]], L5TSR_2001[[1]])

ls_res <- suppressMessages(relnorm(as(L5TSR_1986[[1]], 'SpatialGridDataFrame'),
                                   as(L5TSR_2001[[1]], 'SpatialGridDataFrame'), nperm=0))
ls_res <- raster(ls_res$newimage)

test_that("rastnorm works for RasterLayers", {
          expect_equal(tl_res, ls_res)
})

###############################################################################
# Test on stacks
tl_res_stack <- suppressMessages(normalize(L5TSR_1986, L5TSR_2001))

ls_res_stack <- stack()
for (n in 1:nlayers(L5TSR_2001)) {
    ls_out <- suppressMessages(relnorm(as(L5TSR_1986[[n]], 'SpatialGridDataFrame'),
                                       as(L5TSR_2001[[n]], 'SpatialGridDataFrame'),
                                       nperm=0))
    ls_out <- raster(ls_out$newimage)
    ls_res_stack <- addLayer(ls_res_stack, ls_out)
}

test_that("rastnorm works for RasterStacks", {
          expect_equal(tl_res_stack, ls_res_stack)
})

###############################################################################
# Test on stacks
suppressMessages(library(foreach))
suppressMessages(library(doParallel))

registerDoParallel(2)
tl_res_stack_parallel <- normalize(L5TSR_1986, L5TSR_2001)
stopImplicitCluster()

test_that("rastnorm works in parallel for RasterStacks", {
    expect_equal(tl_res_stack, tl_res_stack_parallel)
})
#
# x <- stack('O:/Data/CI_H_Data/Landsat/BBS/124-063_2000-120_LT5/BBS_124-063_2000-120_L5TSR_tc.envi')
# y <- stack('O:/Data/CI_H_Data/Landsat/BBS/124-063_2000-136_LT5/BBS_124-063_2000-136_L5TSR_tc.envi')
# x_mask <- stack('O:/Data/CI_H_Data/Landsat/BBS/124-063_2000-120_LT5/BBS_124-063_2000-120_L5TSR_masks.envi')
# y_mask <- stack('O:/Data/CI_H_Data/Landsat/BBS/124-063_2000-136_LT5/BBS_124-063_2000-136_L5TSR_masks.envi')
# msk <- (y_mask[[2]] == 255) | (y_mask[[2]] == 2) | (y_mask[[2]] == 4) | 
#        (x_mask[[2]] == 255) | (x_mask[[2]] == 2) | (x_mask[[2]] == 4)
# x_spdf <- as(x[[1]], "SpatialGridDataFrame")
# y_spdf <- as(y[[1]], "SpatialGridDataFrame")
# msk_spdf <- msk
# msk_spdf[msk_spdf == 0] <- NA
# msk_spdf <- as(msk_spdf, "SpatialGridDataFrame")
#
# normed_y_tl <- normalize(x[[1]], y[[2]], msk)
# normed_y_ls <- relnorm(x_spdf, y_spdf, msk_spdf, nperm=0)
#
# normed_y_tl <- normalize(x, y, msk)
#
# normed_y_ls <- raster(normed_y_ls$newimage)
# normed_y_tl
# normed_y_ls
