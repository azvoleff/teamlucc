context("minnart_samp code")

library(landsat)

# Load sample data
DEM_mosaic <- mosaic(ASTER_V002_EAST, ASTER_V002_WEST, fun='mean')
matched_DEM <- match_rasters(L5TSR_1986, DEM_mosaic)
slopeaspect <- slopeasp_seq(matched_DEM)

###############################################################################
# Test that minslope methods match between landsat package 'topocorr' and teamr 
# 'topographic_corr':
sunelev <- 90 - 44.97 # From metadata file
sunazimuth <- 124.37 # From metadata file
L5TSR_1986_b1 <- raster(L5TSR_1986, layer=1)
teamr_tc <- topographic_corr(L5TSR_1986_b1, slopeaspect, sunelev, sunazimuth, 
                             method='minslope')

slope <- as(raster(slopeaspect, layer=1), "SpatialGridDataFrame")
aspect <- as(raster(slopeaspect, layer=2), "SpatialGridDataFrame")
landsat_tc_b1 <- topocorr(as(L5TSR_1986_b1, "SpatialGridDataFrame"), slope, 
                          aspect, sunelev, sunazimuth, method='minslope')
landsat_tc_b1 <- raster(landsat_tc_b1)
names(landsat_tc_b1) <- 'b1tc'

test_that("teamr and landsat minslope match", {
          expect_equal(teamr_tc, expected=landsat_tc_b1)
})

###############################################################################
set.seed(1)

sampleindices <- gridsample(L5TSR_1986_b1, rowmajor=TRUE)
teamr_tc_sample <- topographic_corr(L5TSR_1986_b1, slopeaspect, sunelev, 
                                    sunazimuth, method='minslope', 
                                    sampleindices=sampleindices)

test_that("teamr and landsat minslope match when sampling is used in teamr", {
          expect_equal(teamr_tc_sample, expected=landsat_tc_b1, tolerance=.25)
})

###############################################################################
L5TSR_1986_b2 <- raster(L5TSR_1986, layer=2)
landsat_tc_b2 <- topocorr(as(L5TSR_1986_b2, "SpatialGridDataFrame"), slope, 
                          aspect, sunelev, sunazimuth, method='minslope')
landsat_tc_b2 <- raster(landsat_tc_b2)
names(landsat_tc_b2) <- 'b2tc'
landsat_tc_b1_b2 <- stack(landsat_tc_b1, landsat_tc_b2)

teamr_tc_b1_b2 <- topographic_corr(stack(L5TSR_1986_b1, L5TSR_1986_b2),
                                   slopeaspect, sunelev, sunazimuth, 
                                   method='minslope')

test_that("teamr and landsat minslope match when multiple layers are processed in teamr", {
          expect_equal(teamr_tc_b1_b2, expected=landsat_tc_b1_b2)
})

