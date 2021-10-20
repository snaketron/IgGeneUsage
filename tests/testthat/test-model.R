context("Tests stan model")
cat("Tests stan model \n")

test_that("stan model availability check", {
  
  model.file <- system.file("extdata", "zibb.stan",
                            package = "IgGeneUsage")
  model <- rstan::stan_model(file = model.file,
                             auto_write = TRUE)
  
  expect_equal(object = model@model_name, expected = "zibb")
  expect_is(object = model, class = "stanmodel")
  
  # data(Ig)
  # fit <- DGU(usage.data = Ig,
  #            mcmc.warmup = 250,
  #            mcmc.steps = 500,
  #            mcmc.chains = 1,
  #            mcmc.cores = 1,
  #            hdi.level = 0.95,
  #            adapt.delta = 0.90,
  #            max.treedepth = 10)
  # fit <- fit$glm
  # expect_is(object = fit, class = "stanfit")
})
