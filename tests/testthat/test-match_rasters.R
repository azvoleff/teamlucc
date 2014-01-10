context("match_rasters")

wgs84_string <- "+init=epsg:4326 +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
baseimg <- raster(matrix(1:16, nrow=4), crs=wgs84_string)

###############################################################################
# Test with resampling
matchimg_wgs84 <- raster(matrix(21:36, nrow=4), crs=wgs84_string)
extent(matchimg_wgs84) <- c(0, 1, .5, 1.25)

expected_wgs84 <- raster(matrix(c(23, 27, 31, 35, 24, 28, 32, 36, NA, NA, NA, 
                                  NA, NA, NA, NA, NA), nrow=4, byrow=TRUE), 
                         crs=wgs84_string)
names(expected_wgs84) <- "layer"

test_that("match raster works correctly with resampling", {
    expect_equal(match_rasters(baseimg, matchimg_wgs84, method='ngb'), 
                 expected=expected_wgs84)
})

###############################################################################
# Test with reprojection
matchimg_nad83 <- raster(matrix(21:36, nrow=4), crs="+init=epsg:4629")

expected_nad83 <- raster(matrix(c(21, 25, 29, 33, 22, 26, 30, 34, 23, 27, 31, 
                                  35, 24, 28, 32, 36), nrow=4, byrow=TRUE),
                         crs=wgs84_string)
names(expected_nad83) <- "layer"


get_nad83 <- match_rasters(baseimg, matchimg_nad83, method='ngb')

test_that("match raster works correctly with reprojection", {
    expect_equal(match_rasters(baseimg, matchimg_nad83, method='ngb'), 
                 expected=expected_nad83)
})

