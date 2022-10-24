context("Tests stan model")
cat("Tests stan model \n")


test_that("stan model and DGU availability check", {
  
  expect_equal(object = stanmodels$zibb@model_name, expected = "zibb")
  expect_is(object = stanmodels$zibb, class = "stanmodel")
  
  data(Ig)
  fit <- expect_warning(
    object = DGU(usage.data = Ig,
                 mcmc.warmup = 500,
                 mcmc.steps = 2500,
                 mcmc.chains = 3,
                 mcmc.cores = 1,
                 hdi.level = 0.95,
                 adapt.delta = 0.99,
                 max.treedepth = 10),
    regexp = "cannot compute exact p-value with ties")
  
  fit <- fit$glm
  expect_is(object = fit, class = "stanfit")
})




test_that("stan model and LOO availability check", {
  
  expect_equal(object = stanmodels$zibb@model_name, expected = "zibb")
  expect_is(object = stanmodels$zibb, class = "stanmodel")
  
  data(Ig)
  fit <- expect_error(
    object = LOO(usage.data = Ig,
                 mcmc.warmup = 500,
                 mcmc.steps = 2500,
                 mcmc.chains = 3,
                 mcmc.cores = 1,
                 hdi.level = 0.95,
                 adapt.delta = 0.999,
                 max.treedepth = 10),
    NA)
})
