

diffUsage <- function(usage.data,
                      mcmc.warmup = 500,
                      mcmc.steps = 1500,
                      mcmc.chains = 4,
                      mcmc.cores = 4,
                      hdi.level = 0.95,
                      adapt.delta = 0.99,
                      max.treedepth = 12,
                      dev.model) {


  # check inputs
  checkInput(usage.data = usage.data,
             mcmc.chains = mcmc.chains,
             mcmc.cores = mcmc.cores,
             mcmc.steps = mcmc.steps,
             mcmc.warmup = mcmc.warmup,
             hdi.level = hdi.level)



  usage.data.raw <- usage.data
  # format input usage
  usage.data <- getUsageData(usage = usage.data.raw)


  model <- rstan::stan_model(file = dev.model)
  # model
  # model <- rstan::stan_model(file = "src/stan_files/zib_multiz.stan")


  # stan sampling
  glm <- rstan::sampling(object = model, #stanmodels$model
                         data = usage.data,
                         chains = mcmc.chains,
                         cores = mcmc.cores,
                         iter = mcmc.steps,
                         warmup = mcmc.warmup,
                         control = list(adapt_delta = adapt.delta,
                                        max_treedepth = max.treedepth))


  # get summary
  glm.summary <- summary(object = glm, digits = 4, pars = "beta_gene",
                         prob = c(0.5, (1-hdi.level)/2, 1-(1-hdi.level)/2))
  glm.summary <- glm.summary$summary
  glm.summary <- data.frame(glm.summary)
  colnames(glm.summary) <- c("effect_mean", "effect_mean_se",
                             "effect_sd", "effect_median",
                             "effect_L", "effect_H", "Neff", "Rhat")


  # extract data
  glm.ext <- rstan::extract(object = glm)



  # get pmax
  glm.summary$pmax <- getPmax(glm.ext = glm.ext)



  # Depreceated: bc
  # glm.summary$bc <- getBcStats(glm.ext = glm.ext,
  #                              usage.data = usage.data)



  # add gene id
  glm.summary$gene_name <- usage.data$gene_names



  # ppc
  ppc <- getPpc(glm.ext = glm.ext,
                usage.data = usage.data,
                hdi.level = hdi.level)


  # group ppc
  group.ppc <- getGroupStats(glm.ext = glm.ext,
                             usage.data = usage.data,
                             hdi.level = hdi.level)



  # frequentist tests, merge data
  t.test.stats <- getTTestStats(usage.data = usage.data)
  u.test.stats <- getManUStats(usage.data = usage.data)
  test.stats <- merge(x = t.test.stats, y = u.test.stats,
                      by = "gene_name")


  # result
  result <- list(glm = glm,
                 glm.summary = glm.summary,
                 group.ppc = group.ppc,
                 ppc = ppc,
                 usage.data = usage.data,
                 test.stats = test.stats)

  return (result)
}


