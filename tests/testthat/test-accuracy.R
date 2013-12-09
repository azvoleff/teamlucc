context("accuracy assessment")

# Below test case is from: Olofsson, P., G. M. Foody, S. V. Stehman, and C. E.  
# Woodcock. 2013. Making better use of accuracy data in land change studies: 
# Estimating accuracy and area and quantifying uncertainty using stratified 
# estimation. Remote Sensing of Environment 129:122–131.
ct_olof <- matrix(c(97, 0, 3,
                        3, 279, 18,
                        2, 1, 97), nrow=3, byrow=TRUE)
pop_olof <- c(22353, 1122543, 610228)
pop_ct_olof <- matrix(c(.012, 0, 0.0004,
                        .006, .595, .038,
                        .007, .003, .337), nrow=3, byrow=TRUE)
test_that(".calc_pop_ct works for olof et al. (2013) example", {
          expect_equal(.calc_pop_ct(ct_olof, pop_olof), 
                       expected=pop_ct_olof,
                       tolerance=.003)
})

# This test case is also from Olofsson et al. (2013)
pop_ct_olof_m <- matrix(c(.012, 0, 0.0004, .013, .97,
                          .006, .595, .038, .639, .93,
                          .007, .003, .337, .347, .97,
                          .026, .598, .376, 1, NA,
                          .48, .99, .88, NA, .94), nrow=5, 
                          byrow=TRUE)
dimnames(pop_ct_olof_m) <- list(predicted=c('', '', '', '', 'Producers'),
                                observed=c('', '', '', '', 'Users'))
class(pop_ct_olof_m) <- 'table'
test_that(".add_ct_margins works for Olofsson et al. (2013) example", {
          expect_equal(.add_ct_margins(pop_ct_olof), expected=pop_ct_olof_m, 
                       tolerance=.01)
})

# Below tests are based on "PontiusMatrix26.xlsx" Excel sheet provided by Gil 
# Pontius.
ct_pontius <- matrix(c(11, 12, 13, 14,
                       21, 22, 23, 24,
                       31, 32, 33, 34,
                       41, 42, 43, 44), nrow=4, byrow=TRUE)

test_that(".calc_Q works for Pontius and Millones (2011) example", {
         expect_equal(.calc_Q(.calc_pop_ct(ct_pontius, pop=c(100, 200, 300, 400))), 
                      expected=.182234, tolerance=.00001)
})

test_that(".calc_A works for Pontius and Millones (2011) example", {
         expect_equal(.calc_A(.calc_pop_ct(ct_pontius, pop=c(100, 200, 300, 400))), 
                      expected=.567193, tolerance=.00001)
})
