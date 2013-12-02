library(teamr)
library(testthat)
context("GLCM texture calculations")

set.seed(0)
test_matrix <- matrix(runif(100)*16, nrow=10)
# test_rast <- raster(test_matrix, crs='+init=EPSG:4326')
# test_rast <- cut(test_rast, seq(0, 32))
# writeRaster(test_rast, 'glcm_test_raster.envi', overwrite=TRUE)

textures_envi <- brick('C:/Users/azvoleff/Code/TEAM/teamr-test/glcm_test_raster_ENVI_out_3x3_1-1_32.envi')
textures_envi <- setMinMax(textures_envi)
names(textures_envi) <- c('mean', 'variance', 'homogeneity', 
                          'contrast', 'dissimilarity', 'entropy', 
                          'second_moment', 'correlation')

statistics <- c('mean', 'variance', 'homogeneity', 'contrast', 'dissimilarity', 
                'entropy', 'second_moment', 'correlation')

# Make a function to get 2d matrix from 3d matrix returned by calc_texture_full_image
get_glcm_texture <- function(statistic, ...) {
    if (length(statistic) != 1) {
        stop('length of statistic must be equal to 1')
    }
    text <- calc_texture_full_image(test_matrix, statistic, 32, c(3, 3), c(1, 1))
    return(text[,,1])
}

unload("../../")
install("../../")
library(teamr)

# calc_texture_full_image(test_matrix, 'variance', 32, c(3,3), c(1,1))
# calc_texture_full_image(test_matrix, 'mean', 32, c(3,3), c(1,1))
# calc_texture_full_image(test_matrix, 'homogeneity', 32, c(3,3), c(1,1))
# calc_texture_full_image(test_matrix, 'contrast', 32, c(3,3), c(1,1))

textures_teamr <- stack(apply(calc_texture_full_image(test_matrix, statistics, 
                                                      32, c(3,3), c(1,1)), 3, 
                             raster))
names(textures_teamr) <- statistics

textures_envi
textures_teamr

plot(textures_envi)
plot(textures_teamr)
plot(textures_envi$mean)
plot(textures_teamr$mean)
plot(stack(textures_envi$contrast, textures_teamr$contrast))
plot(stack(textures_envi$variance, textures_teamr$variance))

test_that("GLCM mean is correct", {
    expect_equal(get_glcm_texture('mean'),
                 expected=as.matrix(textures_envi$mean),
                 tolerance=.000001)
})

test_that("GLCM variance is correct", {
    expect_equal(get_glcm_texture('variance'),
                 expected=as.matrix(textures_envi$variance),
                 tolerance=.000001)
})

# test_that("GLCM covariance is correct", {
#     expect_equal(get_glcm_texture('covariance'),
#                  expected=as.matrix(textures_envi$covariance),
#                  tolerance=.000001)
# })

test_that("GLCM homogeneity is correct", {
    expect_equal(get_glcm_texture('homogeneity'),
                 expected=as.matrix(textures_envi$homogeneity),
                 tolerance=.000001)
})

test_that("GLCM contrast is correct", {
    expect_equal(get_glcm_texture('contrast'),
                 expected=as.matrix(textures_envi$contrast),
                 tolerance=.000001)
})

test_that("GLCM dissimilarity is correct", {
    expect_equal(get_glcm_texture('dissimilarity'),
                 expected=as.matrix(textures_envi$dissimilarity),
                 tolerance=.000001)
})

test_that("GLCM entropy is correct", {
    expect_equal(get_glcm_texture('entropy'),
                 expected=as.matrix(textures_envi$entropy),
                 tolerance=.000001)
})

test_that("GLCM second_moment is correct", {
    expect_equal(get_glcm_texture('second_moment'),
                 expected=as.matrix(textures_envi$second_moment),
                 tolerance=.000001)
})

test_that("GLCM correlation is correct", {
    expect_equal(get_glcm_texture('correlation'),
                 expected=as.matrix(textures_envi$correlation),
                 tolerance=.000001)
})
