context("GLCM textures")

suppressMessages(library(raster))

# Make a function to get 2d matrix from 3d matrix returned by glcm
get_teamr_glcm_texture <- function(statistic) {
    if (length(statistic) != 1) {
        stop('length of statistic must be equal to 1')
    }
    texture <- glcm(glcm_test_raster, 32, c(3, 3), c(1, 1), statistic)
    return(getValues(texture))
}

# Test all statistics that are available in EXELIS ENVI match the textures 
# output by teamr
test_that("GLCM mean is correct", {
    expect_equal(get_teamr_glcm_texture('mean_ENVI'),
                 expected=getValues(glcm_test_raster_ENVI_textures$mean),
                 tolerance=.000001)
})

test_that("GLCM variance is correct", {
    expect_equal(get_teamr_glcm_texture('variance_ENVI'),
                 expected=getValues(glcm_test_raster_ENVI_textures$variance),
                 tolerance=.000001)
})

test_that("GLCM homogeneity is correct", {
    expect_equal(get_teamr_glcm_texture('homogeneity'),
                 expected=getValues(glcm_test_raster_ENVI_textures$homogeneity),
                 tolerance=.000001)
})

test_that("GLCM contrast is correct", {
    expect_equal(get_teamr_glcm_texture('contrast'),
                 expected=getValues(glcm_test_raster_ENVI_textures$contrast),
                 tolerance=.000001)
})

test_that("GLCM dissimilarity is correct", {
    expect_equal(get_teamr_glcm_texture('dissimilarity'),
                 expected=getValues(glcm_test_raster_ENVI_textures$dissimilarity),
                 tolerance=.000001)
})

test_that("GLCM entropy is correct", {
    expect_equal(get_teamr_glcm_texture('entropy'),
                 expected=getValues(glcm_test_raster_ENVI_textures$entropy),
                 tolerance=.000001)
})

test_that("GLCM second_moment is correct", {
    expect_equal(get_teamr_glcm_texture('second_moment'),
                 expected=getValues(glcm_test_raster_ENVI_textures$second_moment),
                 tolerance=.000001)
})

test_that("GLCM correlation is correct", {
    expect_equal(get_teamr_glcm_texture('correlation'),
                 expected=getValues(glcm_test_raster_ENVI_textures$correlation),
                 tolerance=.000001)
})

# Test that glcm run on a raster matches the output from running glcm directly 
# on a matrix
glcm_corr_matrix <- glcm(raster::as.matrix(glcm_test_raster), 32, c(3, 3), c(1, 1), 'correlation')
glcm_corr_matrix <- matrix(glcm_corr_matrix, nrow=nrow(glcm_corr_matrix))
test_that("GLCM run on a matrix works correctly", {
    expect_equal(glcm_corr_matrix,
                 expected=raster::as.matrix(glcm_test_raster_ENVI_textures$correlation),
                 tolerance=.000001)
})

glcm_corr_int <- round(glcm(glcm_test_raster, 32, c(3, 3), c(1, 1), 'correlation') * 1000)
test_that("GLCM scaling works correctly when run with scaling and rounding", {
    expect_equal(glcm(glcm_test_raster, 32, c(3, 3), c(1, 1), 'correlation', asinteger=TRUE, scale_factor=1000),
                 expected=glcm_corr_int,
                 tolerance=.000001)
})
