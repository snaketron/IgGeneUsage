context("Tests stan model")
cat("Tests stan model \n")


test_that("stan model and DGU availability check", {
  
  model.file <- system.file("extdata", "zibb.stan",
                            package = "IgGeneUsage")
  model <- rstan::stan_model(file = model.file,
                             auto_write = TRUE)
  
  expect_equal(object = model@model_name, expected = "zibb")
  expect_is(object = model, class = "stanmodel")
  
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
  
  model.file <- system.file("extdata", "zibb.stan",
                            package = "IgGeneUsage")
  model <- rstan::stan_model(file = model.file,
                             auto_write = TRUE)
  
  expect_equal(object = model@model_name, expected = "zibb")
  expect_is(object = model, class = "stanmodel")
  
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
  
  fit <- fit$glm
  expect_is(object = fit, class = "stanfit")
})
