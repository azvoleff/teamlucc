context("train_classifier")

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


################################################################################
# Test training and classifying with models with a factor predictor variable

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
