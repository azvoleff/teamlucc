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

expected_layer <- raster(L5TSR_1986, layer=1) * 10
expected_stack <- stack(raster(L5TSR_1986, layer=1) * 10,
                        raster(L5TSR_1986, layer=2) * 10,
                        raster(L5TSR_1986, layer=3) * 10,
                        raster(L5TSR_1986, layer=4) * 1)

test_that("scaling a RasterLayer works correctly", {
    expect_equal(scale_raster(test_layer), expected=expected_layer)
})

test_that("scaling a RasterStack works correctly", {
    expect_equal(scale_raster(test_stack), expected=expected_stack)
})

test_that("scaling a RasterBrick works correctly", {
    expect_equal(scale_raster(test_brick), expected=expected_stack)
})

test_that("returning scale factors works correctly correctly", {
    expect_equal(scale_raster(test_stack, do_scaling=FALSE), 
                 expected=list(b1=10, b2=10, b3=10, b4=1))
})
