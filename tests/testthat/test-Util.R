context("Tests input rules")
# source("R/Util.R")


test_that("Null input", {
  expect_error(checkInput(usage.data = NULL,
                          mcmc.chains = NULL,
                          mcmc.cores = NULL,
                          mcmc.steps = NULL,
                          mcmc.warmup = NULL,
                          hdi.level = NULL),
               "arguments must be specified")

  expect_error(checkInput(usage.data = NA,
                          mcmc.chains = NA,
                          mcmc.cores = NA,
                          mcmc.steps = NA,
                          mcmc.warmup = NA,
                          hdi.level = NA),
               "arguments must be specified")

  expect_error(checkInput(usage.data = Inf,
                          mcmc.chains = Inf,
                          mcmc.cores = Inf,
                          mcmc.steps = Inf,
                          mcmc.warmup = Inf,
                          hdi.level = Inf),
               "usage.data must be data.frame")

})
