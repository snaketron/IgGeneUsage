context("Tests stan model")
cat("Tests stan model \n")

test_that("stan model and LOO availability check", {
  
  expect_equal(object = stanmodels$dgu@model_name, expected = "dgu")
  expect_equal(object = stanmodels$gu@model_name, expected = "gu")
  expect_is(object = stanmodels$dgu, class = "stanmodel")
  expect_is(object = stanmodels$gu, class = "stanmodel")
})
