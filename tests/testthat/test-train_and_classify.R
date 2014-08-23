context("train_and_classify")

set.seed(0)

train_data <- get_pixels(L5TSR_1986, L5TSR_1986_2001_training, "class_1986", 
                         training=.6)

rf_m <- train_classifier(train_data)
svm_m <- train_classifier(train_data, type="svm")
test_that("can train both SVM and random forest models", {
    expect_equal(class(rf_m$finalModel), "randomForest")
    expect_equal(class(svm_m$finalModel)[[1]], "ksvm")
})

# Test running classify from these models
rf_cl <- classify(L5TSR_1986, rf_m)
svm_cl <- classify(L5TSR_1986, svm_m)
expected_codes <- data.frame(code=c(0, 1),
                             class=c("Forest", "NonForest"))
test_that("can run classify on both SVM and random forest models", {
    expect_equal(names(rf_cl), c("classes", "probs", "codes"))
    expect_equivalent(class(rf_cl$classes), "RasterLayer")
    expect_equivalent(class(rf_cl$probs), "RasterBrick")
    expect_equal(rf_cl$codes, expected_codes)
    expect_equal(names(svm_cl), c("classes", "probs", "codes"))
    expect_equivalent(class(svm_cl$classes), "RasterLayer")
    expect_equivalent(class(svm_cl$probs), "RasterBrick")
    expect_equal(svm_cl$codes, expected_codes)
})

calc_mean_diff <- function(p1, p2) {
    mean(abs(cellStats(p1 - p2, "mean")))
}

test_that("SVM random forest models roughly match", {
    expect_less_than(calc_mean_diff(svm_cl$probs, rf_cl$probs), .05)
    expect_less_than(calc_mean_diff(svm_cl$classes, rf_cl$classes), .07)
})

L5TSR_1986_some_NA <- L5TSR_1986
L5TSR_1986_some_NA$b1[1:100, 1:100] <- NA
rf_cl_some_NA <- classify(L5TSR_1986_some_NA, rf_m)
test_that("classify works with some NAs", {
    expect_equal(cellStats(is.na(rf_cl_some_NA$classes), "sum"), 10015)
})

L5TSR_1986_all_NA <- L5TSR_1986
L5TSR_1986_all_NA$b1 <- NA
rf_cl_all_NA <- classify(L5TSR_1986_all_NA, rf_m)
test_that("classify works with ALL NAs", {
    expect_equal(cellStats(is.na(rf_cl_all_NA$classes), "sum"), 
                 ncell(L5TSR_1986_all_NA))
})

################################################################################
# Test training and classifying with models with a factor predictor variable

set.seed(0)

# Add a "year" layer to L5TSR_1986
L5TSR_1986_year <- L5TSR_1986
L5TSR_1986_year$year <- floor(runif(ncell(L5TSR_1986), 1, 5))
fac_levels <- c(1, 2, 3, 4)

train_data_factor <- get_pixels(L5TSR_1986_year, L5TSR_1986_2001_training, 
                                "class_1986", training=.6)
rf_m_factor <- train_classifier(train_data_factor, 
                                factors=list(year=fac_levels))
svm_m_factor <- train_classifier(train_data_factor, type="svm",
                                 factors=list(year=fac_levels))

# Also check that factor is not encoded unless levels are explicitly supplied 
# to train_classifier
svm_m_factor_not_specified <- train_classifier(train_data_factor, type="svm")

test_that("can train models with factor variables as predictors", {
    expect_equal(class(rf_m$finalModel), "randomForest")
    expect_equal(class(svm_m$finalModel)[[1]], "ksvm")
    expect_equal(svm_m_factor$xlevels, list(year=as.character(fac_levels)))
    expect_equal(rf_m_factor$xlevels, list(year=as.character(fac_levels)))
    expect_equivalent(svm_m_factor_not_specified$xlevels, list())
})

# Test that classify will run on these models:
rf_cl_factor <- classify(L5TSR_1986_year, rf_m_factor, 
                         factors=list(year=fac_levels))
svm_cl_factor <- classify(L5TSR_1986_year, svm_m_factor, 
                          factors=list(year=fac_levels))
test_that("can run classify on both SVM and random forest models", {
    expect_equal(names(rf_cl_factor), c("classes", "probs", "codes"))
    expect_equivalent(class(rf_cl_factor$classes), "RasterLayer")
    expect_equivalent(class(rf_cl_factor$probs), "RasterBrick")
    expect_equal(rf_cl_factor$codes, expected_codes)
    expect_equal(names(svm_cl_factor), c("classes", "probs", "codes"))
    expect_equivalent(class(svm_cl_factor$classes), "RasterLayer")
    expect_equivalent(class(svm_cl_factor$probs), "RasterBrick")
    expect_equal(svm_cl_factor$codes, expected_codes)
})

test_that("SVM random forest models roughly match with factors", {
    expect_less_than(calc_mean_diff(svm_cl$classes, svm_cl_factor$classes), .12)
    expect_less_than(calc_mean_diff(rf_cl$classes, rf_cl_factor$classes), 16)
    expect_less_than(calc_mean_diff(svm_cl$probs, svm_cl_factor$probs), .12)
    expect_less_than(calc_mean_diff(rf_cl$probs, rf_cl_factor$probs), .11)
})
