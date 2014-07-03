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
