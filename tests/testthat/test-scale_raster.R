context("raster_scaling")

test_layer <- raster(L5TSR_1986, layer=1)
test_stack <- stack(L5TSR_1986)
test_brick <- L5TSR_1986

if (!(class(test_layer) == 'RasterLayer')) {
    stop('test layer is not a RasterLayer')
}
if (!(class(test_stack) == 'RasterStack')) {
    stop('test stack is not a RasterStack')
}
if (!(class(test_brick) == 'RasterBrick')) {
    stop('test brick is not a RasterBrick')
}

expected_layer <- raster(L5TSR_1986, layer=1) / 100
expected_layer <- round(expected_layer)
names(expected_layer) <- 'b1'
expected_stack <- L5TSR_1986 / 100
expected_stack <- round(expected_stack)
names(expected_stack) <- c('b1', 'b2', 'b3', 'b4')

test_that("scaling a RasterLayer works correctly", {
    expect_equal(scale_raster(test_layer, max_out=256), expected_layer)
})

test_that("scaling a RasterStack works correctly", {
    expect_equivalent(scale_raster(test_stack, max_out=256), expected_stack)
})

test_that("scaling a RasterBrick works correctly", {
    expect_equivalent(scale_raster(test_brick, max_out=256), expected_stack)
})

test_that("returning scale factors works correctly correctly", {
    expect_equal(scale_raster(test_stack, max_out=256, do_scaling=FALSE), 
                 expected=list(b1=.01, b2=.01, b3=.01, b4=.01))
    expect_equal(scale_raster(test_stack, max_out=32768*2, do_scaling=FALSE), 
                 expected=list(b1=1, b2=1, b3=1, b4=10))
})

suppressMessages(library(foreach))
suppressMessages(library(doParallel))
registerDoParallel(2)
test_that("scaling a RasterStack in parallel works correctly", {
    expect_equivalent(scale_raster(test_stack, max_out=256), expected_stack)
    expect_equal(scale_raster(test_stack, max_out=256, do_scaling=FALSE), 
                 expected=list(b1=.01, b2=.01, b3=.01, b4=.01))
})
stopImplicitCluster()
