context("gridsample")

test_that("gridsample works for different size matrices", {

    testmatrix <- matrix(1:10000, nrow=100)

    gridsample_length <- function(...) {length(gridsample(...))}

    # Test that a single vertical and single horizontal cell can be selected
    expect_that(gridsample_length(testmatrix, horizcells=1, vertcells=1, 
                                  nsamp=10), equals(10))
    # Test that a single horizontal cell can be selected
    expect_that(gridsample_length(testmatrix, horizcells=1, vertcells=10, 
                                  nsamp=10), equals(100))
    # Test that a single vertical cell can be selected
    expect_that(gridsample_length(testmatrix, horizcells=10, vertcells=1, 
                                  nsamp=10), equals(100))
    # Test that a multiple vertical and horizontal cells can be selected
    expect_that(gridsample_length(testmatrix, horizcells=10, vertcells=10, 
                                  nsamp=10), equals(1000))
    # Test that a sample can be drawn over the entire matrix
    expect_that(gridsample_length(testmatrix, horizcells=1, vertcells=1, 
                                  nsamp=10000), equals(10000))
    # Test that an error is thrown if nsamp is too large
    expect_error(gridsample_length(testmatrix, horizcells=10, vertcells=10, 
                                  nsamp=101))
    
    # Test that elements selected by column major and row major indices match 
    set.seed(1)
    colmaj <- gridsample(testmatrix, horizcells=10, vertcells=10, nsamp=10)
    set.seed(1)
    rowmaj <- gridsample(testmatrix, horizcells=10, vertcells=10, nsamp=10, rowmajor=TRUE)
    expect_that(testmatrix[colmaj], equals(t(testmatrix)[rowmaj]))

})
