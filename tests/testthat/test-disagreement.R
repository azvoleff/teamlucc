context("disagreement calculations")

# Test based on "PontiusMatrix26.xlsx" Excel sheet provided by Gil Pontius
ct <- matrix(c(11, 12, 13, 14,
               21, 22, 23, 24,
               31, 32, 33, 34,
               41, 42, 43, 44), nrow=4, byrow=TRUE)

test_that("disagreement is correct for Potinus and Millones (2011) example", {
         expect_equal(disagreement(ct, pop=c(100, 200, 300, 400)), 
                      expected=list(Q=.182234, A=.567193), tolerance=.00001)
})

# Test that calculations work if sample frequencies are used to estimate 
# population frequencies. Note that suppressWarnings is used to suppress the 
# warning message from the disagreement function that sample frequencies are 
# being used to estimate the population frequencies.
test_that("disagreement is correct when assume sample freq is pop freq", {
         expect_equal(suppressWarnings(disagreement(ct)),
                      expected=list(Q=.163636, A=.586364), tolerance=.00001)
})

# Below test case is from: Olofsson, P., G. M. Foody, S. V. Stehman, and C. E.  
# Woodcock. 2013. Making better use of accuracy data in land change studies: 
# Estimating accuracy and area and quantifying uncertainty using stratified 
# estimation. Remote Sensing of Environment 129:122–131.
ct_olofsson <- matrix(c(97, 0, 3,
                        3, 279, 18,
                        2, 1, 7), nrow=3, byrow=TRUE)
pop_olofsson <- c(22353, 1122543, 610228)
test_that("disagreement is correct for Olofsson et al. (2013) example", {
         expect_equal(disagreement(ct_olofsson, pop_olofsson), 
                      expected=list(Q=.075550, A=.073907),
                      tolerance=.00001)
})
