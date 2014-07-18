context("pixel_data class")

###############################################################################
# Test the handling of data by the S4 accuracy methods

set.seed(1)
train_data <- get_pixels(L5TSR_1986, L5TSR_1986_2001_training, "class_1986", 
                         training=.6)

test_that("basic pixel_data methods work", {
    expect_equal(length(train_data), 120)
    expect_equal(levels(train_data), c("Forest", "NonForest"))
    expect_equal(n_train(train_data), 72)
    expect_equal(n_test(train_data), 48)
    expect_equal(n_train(train_data, 'Forest'), 48)
    expect_equal(n_test(train_data, 'Forest'), 20)
})

test_that("pixel_data extract methods works", {
    expect_equal(length(train_data['Forest']), 68)
    expect_equal(length(train_data['NonForest']), 52)
})


test_that("src_name method works", {
    expect_equal(src_name(train_data['Forest']), 
                 paste(with(train_data@pixel_src[train_data@y == 'Forest', ], 
                            paste(src, ID, sep='_'))))
})

all_training <- train_data
training_flag(all_training) <- TRUE
test_that("training_flag methods work for all pixels", {
    expect_equal(training_flag(all_training), rep(TRUE, length(train_data)))
})

training_flag(train_data, 'Forest') <- TRUE
test_that("training_flag methods work for particular class", {
    expect_equal(training_flag(train_data, 'Forest'),
                 rep(TRUE, length(train_data['Forest'])))
})
