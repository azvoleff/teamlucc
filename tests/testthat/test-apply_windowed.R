context('apply_windowed')

x <- raster(L5TSR_1986, layer=1)
min_x <- cellStats(x, 'min')
max_x <- cellStats(x, 'max')

if (require(glcm)) {
    ###############################################################################
    # Test that glcm calculated with apply_windowed matches the output from 
    # running glcm directly
    glcm_apply_windowed <- apply_windowed(x, glcm, edge=c(1, 2), min_x=min_x, 
                                          max_x=max_x,
                                          statistics=c('mean', 'variance'))
    # Force glcm_apply_windowed to be stored in memory so it can compare 
    # identically with glcm_glcm
    glcm_apply_windowed <- glcm_apply_windowed + 0

    glcm_glcm <- glcm(x, statistics=c('mean', 'variance'))
    glcm_glcm <- glcm_glcm + 0

    test_that("glcm calculated with apply_windowed matches glcm calculated 
              directly", {
        expect_equal(glcm_apply_windowed, expected=glcm_glcm, tolerance=1e-7)
    })

    ###############################################################################
    # Test that glcm calculated with apply_windowed matches the output from 
    # running glcm directly when only a single texture measure is calculated

    glcm_apply_windowed_single <- apply_windowed(x, glcm, edge=c(1, 2), 
                                                 min_x=min_x, max_x=max_x, 
                                                 statistics=c('mean'))
    # Force glcm_apply_windowed to be stored in memory so it can compare 
    # identically with glcm_glcm
    glcm_apply_windowed_single <- glcm_apply_windowed_single + 0

    glcm_glcm_single <- glcm(x, statistics=c('mean'))

    test_that("glcm calculated with apply_windowed matches glcm calculated 
              directly", {
        expect_equal(glcm_apply_windowed_single, expected=glcm_glcm_single, 
                     tolerance=1e-7)
    })
}
