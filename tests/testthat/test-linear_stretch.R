context("linear stretch")

test_vals <- seq(1:9)

expected_numeric <- c(0.0000000, 0.1093750, 0.2395833, 0.3697917, 0.5000000, 
                      0.6302083, 0.7604167, 0.8906250, 1.0000000)
expected_matrix <- matrix(expected_numeric, nrow=3)
expected_RasterLayer <- raster(matrix(expected_numeric, nrow=3))
expected_RasterStack <- stack(expected_RasterLayer, expected_RasterLayer)

test_RasterLayer <- raster(matrix(test_vals, nrow=3))
test_RasterStack <- stack(test_RasterLayer, test_RasterLayer + 10)

test_that("stretch works for numeric", {
    expect_equal(linear_stretch(test_vals), expected=expected_numeric, 
                 tolerance=1e-6)
})

test_that("stretch works for matrix", {
    expect_equal(linear_stretch(matrix(test_vals, nrow=3)), 
                 expected=expected_matrix, tolerance=1e-6)
})

test_that("stretch works for RasterLayer", {
    expect_equal(linear_stretch(test_RasterLayer),
                 expected=expected_RasterLayer, tolerance=1e-6)
})

test_that("stretch works for RasterStack", {
    expect_equal(linear_stretch(test_RasterStack),
                 expected=expected_RasterStack, tolerance=1e-6)
})

test_that("stretch works for different max_val", {
    expect_equal(linear_stretch(test_vals, max_val=255),
                 expected=expected_numeric * 255, tolerance=1e-6)
})

expected_numeric_pct10 <- c(0.00000, 0.03125, 0.18750, 0.34375, 0.50000, 0.65625, 
                      0.81250, 0.96875, 1.00000)
test_that("stretch works for different pct", {
    expect_equal(linear_stretch(test_vals, pct=10),
                 expected=expected_numeric_pct10, tolerance=1e-6)
})

test_that("stretch works for different pct and max_val", {
    expect_equal(linear_stretch(test_vals, pct=10, max_val=255),
                 expected=expected_numeric_pct10 * 255, tolerance=1e-6)
})
