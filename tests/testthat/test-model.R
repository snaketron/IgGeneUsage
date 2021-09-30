context("Tests stan model")

test_that("stan model availability check", {
  
  model.file <- system.file("extdata", "zibb.stan",
                            package = "IgGeneUsage")
  model <- rstan::stan_model(file = model.file,
                             auto_write = TRUE)
  
  expect_equal(object = model@model_name, expected = "zibb")
  expect_is(object = model, class = "stanmodel")
  fit <- rstan::sampling(object = model)
  expect_is(object = fit, class = "stanfit")
  
  data(Ig)
  fit <- DGU(usage.data = Ig,
             mcmc.warmup = 500,
             mcmc.steps = 1500,
             mcmc.chains = 1,
             mcmc.cores = 1,
             hdi.level = 0.95,
             adapt.delta = 0.95,
             max.treedepth = 13)
  fit <- fit$glm
  expect_is(object = fit, class = "stanfit")
})
