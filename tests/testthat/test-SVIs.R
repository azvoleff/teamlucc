context("SVI equation checks")

test_that("NDVI equation is correct", {
    expect_that(NDVI(red=10, nir=6), equals(-.25))
    expect_that(NDVI(red=3, nir=7), equals(.4))
})

test_that("EVI equation is correct", {
    expect_that(EVI(blue=1, red=5, nir=14), equals(.6))
    expect_that(EVI(blue=13, red=5, nir=4), equals(.04))
})

test_that("MSAVI2 equation is correct", {
    expect_that(MSAVI2(red=7, nir=3), equals(-1))
    expect_that(MSAVI2(red=18, nir=5), equals(-2))
})

test_that("ARVI equation is correct", {
    expect_that(ARVI(blue=1, red=5, nir=16), equals(.2))
    expect_that(ARVI(blue=16, red=3, nir=4), equals(3))
})
