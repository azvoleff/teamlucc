context("topographic_corr")

suppressMessages(library(landsat))

# Load sample data
L5TSR_1986_b1 <- raster(L5TSR_1986, layer=1)
L5TSR_1986_b2 <- raster(L5TSR_1986, layer=2)
DEM_mosaic <- mosaic(ASTER_V002_EAST, ASTER_V002_WEST, fun='mean')
matched_DEM <- match_rasters(L5TSR_1986, DEM_mosaic)
slopeaspect <- terrain(matched_DEM, opt=c('slope', 'aspect'))
slope_deg <- as(slopeaspect$slope * (180/pi), "SpatialGridDataFrame")
aspect_deg  <- as(slopeaspect$aspect * (180/pi), "SpatialGridDataFrame")
sunelev <- 90 - 44.97 # From metadata file
sunazimuth <- 124.37 # From metadata file

###############################################################################
# Test that minslope methods match between landsat package 'topocorr' and teamr 
# 'topographic_corr':
teamr_tc_b1 <- topographic_corr(L5TSR_1986_b1, slopeaspect, sunelev, sunazimuth, 
                                method='minslope')

landsat_tc_b1 <- topocorr(as(L5TSR_1986_b1, "SpatialGridDataFrame"), slope_deg,
                          aspect_deg, sunelev, sunazimuth, method='minslope')
landsat_tc_b1 <- raster(landsat_tc_b1)
names(landsat_tc_b1) <- 'b1tc'

test_that("teamr and landsat minslope match", {
          expect_equal(teamr_tc_b1, expected=landsat_tc_b1)
})

###############################################################################
set.seed(1)
sampleindices <- gridsample(L5TSR_1986_b1, rowmajor=TRUE)
teamr_tc_b1_sample <- topographic_corr(L5TSR_1986_b1, slopeaspect, sunelev, 
                                       sunazimuth, method='minslope', 
                                       sampleindices=sampleindices)

test_that("teamr and landsat minslope match when sampling is used in teamr", {
          expect_equal(teamr_tc_b1_sample, expected=landsat_tc_b1, tolerance=.25)
})

###############################################################################
landsat_tc_b2 <- topocorr(as(L5TSR_1986_b2, "SpatialGridDataFrame"), slope_deg, 
                          aspect_deg, sunelev, sunazimuth, method='minslope')
landsat_tc_b2 <- raster(landsat_tc_b2)
names(landsat_tc_b2) <- 'b2tc'
landsat_tc_b1_b2 <- stack(landsat_tc_b1, landsat_tc_b2)

teamr_tc_b1_b2 <- topographic_corr(stack(L5TSR_1986_b1, L5TSR_1986_b2),
                                   slopeaspect, sunelev, sunazimuth, 
                                   method='minslope')

test_that("teamr and landsat minslope match when multiple layers are processed in teamr", {
          expect_equal(teamr_tc_b1_b2, expected=landsat_tc_b1_b2)
})

###############################################################################
# Test that minnaert_full methods match between landsat package 'topocorr' and 
# teamr 'topographic_corr' when using full image.
teamr_minnaert <- topographic_corr(L5TSR_1986_b1, slopeaspect, sunelev, 
                                   sunazimuth, method='minnaert_full', 
                                   slopeclass=c(1, 5, 10, 15, 20, 25, 30, 45))

landsat_minnaert <- minnaert(as(L5TSR_1986_b1, "SpatialGridDataFrame"), slope_deg, 
                             aspect_deg, sunelev, sunazimuth)
landsat_minnaert <- raster(landsat_minnaert$minnaert)
names(landsat_minnaert) <- 'b1tc'

# Need a tolerance since teamr uses 'bam' instead of 'gam'
test_that("teamr minnaert and landsat minnaert match", {
          expect_equal(teamr_minnaert, expected=landsat_minnaert, tolerance=.1)
})

###############################################################################
# Test that minnaert_full methods match between landsat package 'topocorr' and 
# teamr 'topographic_corr' when using resampling.
set.seed(0)
sampleindices <- gridsample(L5TSR_1986_b1, rowmajor=TRUE)
teamr_minnaert_sample <- topographic_corr(L5TSR_1986_b1, slopeaspect, sunelev, 
                                          sunazimuth, method='minnaert_full',
                                          sampleindices=sampleindices,
                                          slopeclass=c(1, 5, 10, 15, 20, 25, 30, 45))

test_that("teamr minnaert sample and landsat minnaert match", {
          expect_equal(teamr_minnaert_sample, expected=landsat_minnaert, 
                       tolerance=.25)
})

###############################################################################
# Test that minnaert_full methods match when teamr minnaert_full runs 
# sequentially or in parallel.

## Commented out as R CMD CHECK hangs when running parallel tests.
# suppressMessages(library(spatial.tools))
# set.seed(0)
# sampleindices <- gridsample(L5TSR_1986_b1, rowmajor=TRUE)
# sfQuickInit(2)
#
# teamr_minnaert_sample_b1b2 <- topographic_corr(stack(L5TSR_1986_b1, L5TSR_1986_b2),
#                                               slopeaspect, sunelev, sunazimuth, 
#                                               method='minnaert_full',
#                                               sampleindices=sampleindices)
#
# teamr_minnaert_sample_b1b2_par <- topographic_corr(stack(L5TSR_1986_b1, L5TSR_1986_b2),
#                                               slopeaspect, sunelev, sunazimuth, 
#                                               method='minnaert_full',
#                                               sampleindices=sampleindices,
#                                               inparallel=TRUE)
# sfQuickStop(2)
#
# test_that("teamr minnaert sample and landsat minnaert match", {
#           expect_equal(teamr_minnaert_sample_b1b2, expected=teamr_minnaert_sample_b1b2_par)
# })
