context("disagreement calculations")

# Test based on "PontiusMatrix26.xlsx" Excel sheet provided by Gil Pontius
ct <- matrix(c(11, 12, 13, 14,
               21, 22, 23, 24,
               31, 32, 33, 34,
               41, 42, 43, 44), nrow=4, byrow=TRUE)

test_that("disagreement is correct for ", {
         expect_equal(disagreement(ct, pop=c(100, 200, 300, 400)), 
                      expected=list(Q=18, A=57))
})

test_that("disagreement is correct for ", {
         expect_equal(disagreement(ct, pop=c(100, 200, 300, 400)), 
                      expected=list(Q=16, A=59))
})
