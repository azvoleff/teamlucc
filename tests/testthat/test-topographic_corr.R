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
# cosine
tl_tc_b1_cosine <- topographic_corr(L5TSR_1986_b1, slopeaspect, sunelev, 
                             sunazimuth, method='cosine')

ls_tc_b1_cosine <- topocorr(as(L5TSR_1986_b1, "SpatialGridDataFrame"), slope_deg,
                     aspect_deg, sunelev, sunazimuth, method='cosine')
ls_tc_b1_cosine <- raster(ls_tc_b1_cosine)
names(ls_tc_b1_cosine) <- 'b1tc'

test_that("teamlucc and landsat cosine match", {
    expect_equal(tl_tc_b1_cosine, expected=ls_tc_b1_cosine)
})

###############################################################################
# improvedcosine
tl_tc_b1_improvedcosine <- topographic_corr(L5TSR_1986_b1, slopeaspect, sunelev, 
                             sunazimuth, method='improvedcosine')

ls_tc_b1_improvedcosine <- topocorr(as(L5TSR_1986_b1, "SpatialGridDataFrame"), slope_deg,
                     aspect_deg, sunelev, sunazimuth, method='improvedcosine')
ls_tc_b1_improvedcosine <- raster(ls_tc_b1_improvedcosine)
names(ls_tc_b1_improvedcosine) <- 'b1tc'

test_that("teamlucc and landsat improvedcosine match", {
    expect_equal(tl_tc_b1_improvedcosine, expected=ls_tc_b1_improvedcosine)
})

###############################################################################
# minnaert
tl_tc_b1_minnaert <- topographic_corr(L5TSR_1986_b1, slopeaspect, sunelev, 
                             sunazimuth, method='minnaert')

ls_tc_b1_minnaert <- topocorr(as(L5TSR_1986_b1, "SpatialGridDataFrame"), slope_deg,
                     aspect_deg, sunelev, sunazimuth, method='minnaert')
ls_tc_b1_minnaert <- raster(ls_tc_b1_minnaert)
names(ls_tc_b1_minnaert) <- 'b1tc'

test_that("teamlucc and landsat minnaert match", {
    expect_equal(tl_tc_b1_minnaert, expected=ls_tc_b1_minnaert)
})

###############################################################################
# minslope
tl_tc_b1_minslope <- topographic_corr(L5TSR_1986_b1, slopeaspect, sunelev, 
                                      sunazimuth, method='minslope')

ls_tc_b1_minslope <- topocorr(as(L5TSR_1986_b1, "SpatialGridDataFrame"), 
                              slope_deg, aspect_deg, sunelev, sunazimuth, 
                              method='minslope')
ls_tc_b1_minslope <- raster(ls_tc_b1_minslope)
names(ls_tc_b1_minslope) <- 'b1tc'

test_that("teamlucc and landsat minslope match", {
    expect_equal(tl_tc_b1_minslope, expected=ls_tc_b1_minslope)
})

###############################################################################
# ccorrection
tl_tc_b1_ccorrection <- topographic_corr(L5TSR_1986_b1, slopeaspect, sunelev, 
                                         sunazimuth, method='ccorrection')

ls_tc_b1_ccorrection <- topocorr(as(L5TSR_1986_b1, "SpatialGridDataFrame"), 
                                 slope_deg, aspect_deg, sunelev, sunazimuth, 
                                 method='ccorrection')
ls_tc_b1_ccorrection <- raster(ls_tc_b1_ccorrection)
names(ls_tc_b1_ccorrection) <- 'b1tc'

test_that("teamlucc and landsat ccorrection match", {
    expect_equal(tl_tc_b1_ccorrection, expected=ls_tc_b1_ccorrection)
})

###############################################################################
# gamma
tl_tc_b1_gamma <- topographic_corr(L5TSR_1986_b1, slopeaspect, sunelev, 
                                   sunazimuth, method='gamma')

ls_tc_b1_gamma <- topocorr(as(L5TSR_1986_b1, "SpatialGridDataFrame"), 
                           slope_deg, aspect_deg, sunelev, sunazimuth, 
                           method='gamma')
ls_tc_b1_gamma <- raster(ls_tc_b1_gamma)
names(ls_tc_b1_gamma) <- 'b1tc'

test_that("teamlucc and landsat gamma match", {
    expect_equal(tl_tc_b1_gamma, expected=ls_tc_b1_gamma)
})

###############################################################################
# SCS
tl_tc_b1_SCS <- topographic_corr(L5TSR_1986_b1, slopeaspect, sunelev, 
                                 sunazimuth, method='SCS')

ls_tc_b1_SCS <- topocorr(as(L5TSR_1986_b1, "SpatialGridDataFrame"), slope_deg, 
                         aspect_deg, sunelev, sunazimuth, method='SCS')
ls_tc_b1_SCS <- raster(ls_tc_b1_SCS)
names(ls_tc_b1_SCS) <- 'b1tc'

test_that("teamlucc and landsat SCS match", {
    expect_equal(tl_tc_b1_SCS, expected=ls_tc_b1_SCS)
})

###############################################################################
# illumination
tl_tc_b1_illumination <- topographic_corr(L5TSR_1986_b1, slopeaspect, sunelev, 
                                          sunazimuth, method='illumination')

ls_tc_b1_illumination <- topocorr(as(L5TSR_1986_b1, "SpatialGridDataFrame"), 
                                  slope_deg, aspect_deg, sunelev, sunazimuth, 
                                  method='illumination')
ls_tc_b1_illumination <- raster(ls_tc_b1_illumination)
names(ls_tc_b1_illumination) <- 'b1tc'

test_that("teamlucc and landsat illumination match", {
    expect_equal(tl_tc_b1_illumination, expected=ls_tc_b1_illumination)
})

###############################################################################
# minnaert with sampling
set.seed(1)
sampleindices <- gridsample(L5TSR_1986_b1, rowmajor=TRUE)
tl_tc_b1_minnaert_sample <- topographic_corr(L5TSR_1986_b1, slopeaspect, 
                                             sunelev, sunazimuth, 
                                             method='minnaert', 
                                             sampleindices=sampleindices)

test_that("teamlucc and landsat minnaert match when sampling is used in teamlucc", {
    expect_equal(tl_tc_b1_minnaert_sample, expected=ls_tc_b1_minnaert, 
                 tolerance=.25)
})

###############################################################################
# minslope with sampling
set.seed(1)
sampleindices <- gridsample(L5TSR_1986_b1, rowmajor=TRUE)
tl_tc_b1_minslope_sample <- topographic_corr(L5TSR_1986_b1, slopeaspect, 
                                             sunelev, sunazimuth, 
                                             method='minslope', 
                                             sampleindices=sampleindices)

test_that("teamlucc and landsat minslope match when sampling is used in teamlucc", {
    expect_equal(tl_tc_b1_minslope_sample, expected=ls_tc_b1_minslope, 
                 tolerance=.25)
})

###############################################################################
# ccorrection with sampling
set.seed(1)
sampleindices <- gridsample(L5TSR_1986_b1, rowmajor=TRUE)
tl_tc_b1_ccorrection_sample <- topographic_corr(L5TSR_1986_b1, slopeaspect, 
                                                sunelev, sunazimuth, 
                                                method='ccorrection', 
                                                sampleindices=sampleindices)

test_that("teamlucc and landsat ccorrection match when sampling is used in teamlucc", {
    expect_equal(tl_tc_b1_ccorrection_sample, expected=ls_tc_b1_ccorrection, 
                 tolerance=.25)
})

###############################################################################
# minslope on raster stack
ls_tc_b2_minslope <- topocorr(as(L5TSR_1986_b2, "SpatialGridDataFrame"), 
                              slope_deg, aspect_deg, sunelev, sunazimuth, 
                              method='minslope')
ls_tc_b2_minslope <- raster(ls_tc_b2_minslope)
names(ls_tc_b2_minslope) <- 'b2tc'
ls_tc_b1b2_minslope <- stack(ls_tc_b1_minslope, ls_tc_b2_minslope)

tl_tc_b1b2_minslope <- topographic_corr(stack(L5TSR_1986_b1, L5TSR_1986_b2), 
                                        slopeaspect, sunelev, sunazimuth, 
                                        method='minslope')

test_that("teamlucc and landsat minslope match when multiple layers are processed in teamlucc", {
    expect_equal(tl_tc_b1b2_minslope, expected=ls_tc_b1b2_minslope)
})

###############################################################################
# minnaert_full
tl_min <- topographic_corr(L5TSR_1986_b1, slopeaspect, sunelev, 
                           sunazimuth, method='minnaert_full')

ls_min <- minnaert(as(L5TSR_1986_b1, "SpatialGridDataFrame"), slope_deg, 
                   aspect_deg, sunelev, sunazimuth)
ls_min <- raster(ls_min$minnaert)
names(ls_min) <- 'b1tc'

# Need a tolerance since teamlucc uses 'bam' instead of 'gam'
test_that("teamlucc minnaert and landsat minnaert match", {
    expect_equal(tl_min, expected=ls_min, tolerance=.1)
})

###############################################################################
# minnaert_full with sampling
set.seed(1)
sampleindices <- gridsample(L5TSR_1986_b1, rowmajor=TRUE)
tl_min_sample <- topographic_corr(L5TSR_1986_b1, slopeaspect, sunelev, 
                                  sunazimuth, method='minnaert_full',
                                  sampleindices=sampleindices)

test_that("teamlucc minnaert sample and landsat minnaert match", {
    expect_equal(tl_min_sample, expected=ls_min, tolerance=.25)
})

###############################################################################
# minnaert_full in parallel
tl_min_b1b2_seq <- topographic_corr(stack(L5TSR_1986_b1, L5TSR_1986_b2),
                                    slopeaspect, sunelev, sunazimuth, 
                                    method='minnaert_full')

suppressMessages(library(foreach))
suppressMessages(library(doParallel))
registerDoParallel(2)
tl_min_b1b2_par <- topographic_corr(stack(L5TSR_1986_b1, L5TSR_1986_b2),
                                    slopeaspect, sunelev, sunazimuth, 
                                    method='minnaert_full')
stopImplicitCluster()

test_that("teamlucc minnaert sample and landsat minnaert match when run in parallel", {
    expect_equal(tl_min_b1b2_seq, expected=tl_min_b1b2_par)
})

###############################################################################
# minnaert_full in parallel with sampling
set.seed(1)
sampleindices <- gridsample(L5TSR_1986_b1, rowmajor=TRUE)
tl_min_sample_b1b2_seq <- topographic_corr(stack(L5TSR_1986_b1, L5TSR_1986_b2),
                                           slopeaspect, sunelev, sunazimuth, 
                                           method='minnaert_full',
                                           sampleindices=sampleindices)

suppressMessages(library(foreach))
suppressMessages(library(doParallel))
registerDoParallel(2)
tl_min_sample_b1b2_par <- topographic_corr(stack(L5TSR_1986_b1, L5TSR_1986_b2),
                                           slopeaspect, sunelev, sunazimuth, 
                                           method='minnaert_full',
                                           sampleindices=sampleindices)
stopImplicitCluster()
