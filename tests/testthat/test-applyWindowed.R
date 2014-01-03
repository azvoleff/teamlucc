context('applyWindowed')

x <- raster(L5TSR_1986, layer=1)
min_x <- cellStats(x, 'min')
max_x <- cellStats(x, 'max')

###############################################################################
# Test that glcm calculated with applyWindowed matches the output from running 
# glcm directly
#TODO: below line is temporary
load_all('../')
glcm_applyWindowed <- applyWindowed(x, glcm, edge=c(1, 2), min_x=min_x, 
                                    max_x=max_x)
# Force glcm_applyWindowed to be stored in memory so it can compare 
# identically with glcm_glcm
glcm_applyWindowed <- glcm_applyWindowed + 0

glcm_glcm <- glcm(x)

plot(stack(glcm_glcm$glcm_correlation, glcm_applyWindowed$glcm_correlation))

test_that("glcm calculated with applyWindowed matches glcm calculated directly", {
    expect_equal(glcm_applyWindowed, expected=glcm_glcm)
})

###############################################################################
# Test that glcm calculated with applyWindowed matches the output from running 
# glcm directly when only a single texture measure is calculated

glcm_applyWindowed_single <- applyWindowed(x, glcm, edge=c(1, 2), min_x=min_x, 
                                           max_x=max_x, statistics=c('mean'))
# Force glcm_applyWindowed to be stored in memory so it can compare 
# identically with glcm_glcm
glcm_applyWindowed_single <- glcm_applyWindowed_single + 0

glcm_glcm_single <- glcm(x)

test_that("glcm calculated with applyWindowed matches glcm calculated directly", {
    expect_equal(glcm_applyWindowed_single, expected=glcm_glcm_single)
})

###############################################################################
# Test that glcm calculated with applyWindowed matches the output from running 
# glcm directly when no edge is used
