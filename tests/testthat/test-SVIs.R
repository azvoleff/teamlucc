context("SVIs - NDVI")

NDVI_red_numeric <- c(10, 3)
NDVI_nir_numeric <- c(6, 7)
NDVI_expected_numeric <- c(-.25, .4)
test_that("NDVI works for numeric", {
    expect_that(NDVI(red=NDVI_red_numeric, nir=NDVI_nir_numeric), 
                equals(NDVI_expected_numeric))
})

NDVI_red_matrix <- matrix(NDVI_red_numeric)
NDVI_nir_matrix <- matrix(NDVI_nir_numeric)
NDVI_expected_matrix <- matrix(NDVI_expected_numeric)
test_that("NDVI from matrix works", {
    expect_that(NDVI(red=NDVI_red_matrix, nir=NDVI_nir_matrix), 
                equals(NDVI_expected_matrix))
})

NDVI_red_rast <- raster(NDVI_red_matrix)
NDVI_nir_rast <- raster(NDVI_nir_matrix)
NDVI_expected_rast <- raster(NDVI_expected_matrix)
test_that("NDVI from raster works", {
    expect_that(NDVI(red=NDVI_red_rast, nir=NDVI_nir_rast), 
                equals(NDVI_expected_rast))
})

###############################################################################
context("SVIs - EVI")

EVI_blue_numeric <- c(1, 13)
EVI_red_numeric <- c(5, 5)
EVI_nir_numeric <- c(14, 4)
EVI_expected_numeric <- c(.6, .04)
test_that("EVI works for numeric", {
    expect_that(EVI(blue=EVI_blue_numeric, red=EVI_red_numeric, 
                    nir=EVI_nir_numeric), equals(EVI_expected_numeric))
})

EVI_blue_matrix <- matrix(EVI_blue_numeric)
EVI_red_matrix <- matrix(EVI_red_numeric)
EVI_nir_matrix <- matrix(EVI_nir_numeric)
EVI_expected_matrix <- matrix(EVI_expected_numeric)
test_that("EVI from matrix works", {
    expect_that(EVI(blue=EVI_blue_matrix, red=EVI_red_matrix, 
                    nir=EVI_nir_matrix), equals(EVI_expected_matrix))
})

EVI_blue_rast <- raster(EVI_blue_matrix)
EVI_red_rast <- raster(EVI_red_matrix)
EVI_nir_rast <- raster(EVI_nir_matrix)
EVI_expected_rast <- raster(EVI_expected_matrix)
test_that("EVI from raster works", {
    expect_that(EVI(blue=EVI_blue_rast, red=EVI_red_rast, nir=EVI_nir_rast), 
                equals(EVI_expected_rast))
})


###############################################################################
context("SVIs - MSAVI2")

MSAVI2_red_numeric <- c(7, 18)
MSAVI2_nir_numeric <- c(3, 5)
MSAVI2_expected_numeric <- c(-1, -2)
test_that("MSAVI2 works for numeric", {
    expect_that(MSAVI2(red=MSAVI2_red_numeric, nir=MSAVI2_nir_numeric), 
                equals(MSAVI2_expected_numeric))
})

MSAVI2_red_matrix <- matrix(MSAVI2_red_numeric)
MSAVI2_nir_matrix <- matrix(MSAVI2_nir_numeric)
MSAVI2_expected_matrix <- matrix(MSAVI2_expected_numeric)
test_that("MSAVI2 from matrix works", {
    expect_that(MSAVI2(red=MSAVI2_red_matrix, nir=MSAVI2_nir_matrix), 
                equals(MSAVI2_expected_matrix))
})

MSAVI2_red_rast <- raster(MSAVI2_red_matrix)
MSAVI2_nir_rast <- raster(MSAVI2_nir_matrix)
MSAVI2_expected_rast <- raster(MSAVI2_expected_matrix)
test_that("MSAVI2 from raster works", {
    expect_that(MSAVI2(red=MSAVI2_red_rast, nir=MSAVI2_nir_rast), 
                equals(MSAVI2_expected_rast))
})

###############################################################################
context("SVIs - ARVI")

ARVI_blue_numeric <- c(1, 16)
ARVI_red_numeric <- c(5, 3)
ARVI_nir_numeric <- c(16, 4)
ARVI_expected_numeric <- c(.2, 3)
test_that("ARVI works for numeric", {
    expect_that(ARVI(blue=ARVI_blue_numeric, red=ARVI_red_numeric, 
                     nir=ARVI_nir_numeric), equals(ARVI_expected_numeric))
})

ARVI_blue_matrix <- matrix(ARVI_blue_numeric)
ARVI_red_matrix <- matrix(ARVI_red_numeric)
ARVI_nir_matrix <- matrix(ARVI_nir_numeric)
ARVI_expected_matrix <- matrix(ARVI_expected_numeric)
test_that("ARVI from matrix works", {
    expect_that(ARVI(blue=ARVI_blue_matrix, red=ARVI_red_matrix, 
                     nir=ARVI_nir_matrix), equals(ARVI_expected_matrix))
})

ARVI_blue_rast <- raster(ARVI_blue_matrix)
ARVI_red_rast <- raster(ARVI_red_matrix)
ARVI_nir_rast <- raster(ARVI_nir_matrix)
ARVI_expected_rast <- raster(ARVI_expected_matrix)
test_that("ARVI from raster works", {
    expect_that(ARVI(blue=ARVI_blue_rast, red=ARVI_red_rast, 
                     nir=ARVI_nir_rast), equals(ARVI_expected_rast))
})

