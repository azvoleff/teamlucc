context("GLCM texture calculations")

# Make a function to get 2d matrix from 3d matrix returned by glcm
get_glcm_texture <- function(statistic, ...) {
    if (length(statistic) != 1) {
        stop('length of statistic must be equal to 1')
    }
    texture <- glcm(glcm_test_raster, 32, c(3, 3), c(1, 1), statistic, ...)
    return(as.array(texture)[,,1])
}

# Test all statistics that are available in EXELIS ENVI match the textures 
# output by teamr
statistics <- c('mean_ENVI', 'variance_ENVI', 'homogeneity', 'contrast', 
                'dissimilarity', 'entropy', 'second_moment', 'correlation')

# glcm_test_raster_ENVI_textures contains the results from running ENVI on the 
# glcm_test_raster. The below tests ensure that the results from teamr match 
# those from EXELIS ENVI.
textures_teamr <- glcm(glcm_test_raster, 32, c(3, 3), c(1, 1), statistics)

test_that("GLCM mean is correct", {
    expect_equal(get_glcm_texture('mean_ENVI'),
                 expected=as.matrix(glcm_test_raster_ENVI_textures$mean),
                 tolerance=.000001)
})

test_that("GLCM variance is correct", {
    expect_equal(get_glcm_texture('variance_ENVI'),
                 expected=as.matrix(glcm_test_raster_ENVI_textures$variance),
                 tolerance=.000001)
})

test_that("GLCM homogeneity is correct", {
    expect_equal(get_glcm_texture('homogeneity'),
                 expected=as.matrix(glcm_test_raster_ENVI_textures$homogeneity),
                 tolerance=.000001)
})

test_that("GLCM contrast is correct", {
    expect_equal(get_glcm_texture('contrast'),
                 expected=as.matrix(glcm_test_raster_ENVI_textures$contrast),
                 tolerance=.000001)
})

test_that("GLCM dissimilarity is correct", {
    expect_equal(get_glcm_texture('dissimilarity'),
                 expected=as.matrix(glcm_test_raster_ENVI_textures$dissimilarity),
                 tolerance=.000001)
})

test_that("GLCM entropy is correct", {
    expect_equal(get_glcm_texture('entropy'),
                 expected=as.matrix(glcm_test_raster_ENVI_textures$entropy),
                 tolerance=.000001)
})

test_that("GLCM second_moment is correct", {
    expect_equal(get_glcm_texture('second_moment'),
                 expected=as.matrix(glcm_test_raster_ENVI_textures$second_moment),
                 tolerance=.000001)
})

test_that("GLCM correlation is correct", {
    expect_equal(get_glcm_texture('correlation'),
                 expected=as.matrix(glcm_test_raster_ENVI_textures$correlation),
                 tolerance=.000001)
})
