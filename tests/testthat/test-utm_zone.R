context('utm_zone')

test_that("works for UTM north zone in western hemisphere", {
    expect_equal(utm_zone(-45, 10), "23N")
})

test_that("works for UTM north zone proj4string in western hemisphere", {
    expect_equal(utm_zone(-45, 10, proj4string=TRUE), "+init=epsg:32623")
})

test_that("works for UTM south zone in western hemisphere", {
    expect_equal(utm_zone(-45, -10), "23S")
})

test_that("works for UTM south zone proj4string in western hemisphere", {
    expect_equal(utm_zone(-45, -10, proj4string=TRUE), "+init=epsg:32723")
})

test_that("works for UTM zone including 0 latitude", {
    expect_equal(utm_zone(-45, 0), "23N")
})

test_that("works for UTM zone proj4string including 0 latitude", {
    expect_equal(utm_zone(-45, 0, proj4string=TRUE), "+init=epsg:32623")
})

test_that("works for UTM zone in eastern hemisphere", {
    expect_equal(utm_zone(23, 10), "34N")
})

test_that("works for UTM zone proj4string in eastern hemisphere", {
    expect_equal(utm_zone(23, 10, proj4string=TRUE), "+init=epsg:32634")
})

