\docType{data}
\name{glcm_test_raster}
\alias{glcm_test_raster}
\title{Randomly generated 100x100 test image}
\description{
Used in testing the output from the GLCM texture statistics
C++ code.
}
\examples{
# The image was generated with the following code:
set.seed(0)
test_matrix <- matrix(runif(100)*32, nrow=10)
glcm_test_raster <- raster(test_matrix, crs='+init=EPSG:4326')
glcm_test_raster <- cut(glcm_test_raster, seq(0, 32))
}

