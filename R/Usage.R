

diffUsage <- function(usage.data,
                      mcmc.warmup = 500,
                      mcmc.steps = 1500,
                      mcmc.chains = 4,
                      mcmc.cores = 4,
                      hdi.level = 0.95,
                      adapt.delta = 0.99,
                      max.treedepth = 12,
                      dev.model) {


  # TODO check ref.group


  # check inputs
  checkInput(usage.data = usage.data,
             mcmc.chains = mcmc.chains,
             mcmc.cores = mcmc.cores,
             mcmc.steps = mcmc.steps,
             mcmc.warmup = mcmc.warmup,
             hdi.level = hdi.level)



  # format input usage
  usage.data.raw <- usage.data
  usage.data <- getUsageData(usage = usage.data.raw)


  # contrast
  contrast <- paste("Contrast: ", unique(usage.data$Xorg[usage.data$X == 1]),
                    " - ", unique(usage.data$Xorg[usage.data$X == -1]),
                    sep = '')


  # model
  model <- rstan::stan_model(file = dev.model)
  # model <- stanmodels$zibb_model




  # setup control list
  control.list <- list(adapt_delta = adapt.delta,
                       max_treedepth = max.treedepth)


  # stan sampling
  glm <- rstan::sampling(object = model,
                         data = usage.data,
                         chains = mcmc.chains,
                         cores = mcmc.cores,
                         iter = mcmc.steps,
                         warmup = mcmc.warmup,
                         refresh = 500, # let user define as well?
                         control = control.list)

  # get summary
  glm.summary <- rstan::summary(object = glm, digits = 4,
                                pars = "beta_gene",
                         prob = c(0.5, (1-hdi.level)/2,
                                  1-(1-hdi.level)/2))
  glm.summary <- glm.summary$summary
  glm.summary <- data.frame(glm.summary)
  colnames(glm.summary) <- c("effect_mean", "effect_mean_se",
                             "effect_sd", "effect_median",
                             "effect_L", "effect_H",
                             "Neff", "Rhat")
  glm.summary[, c("Rhat", "Neff")] <- NULL
  glm.summary$contrast <- contrast


  # extract data
  glm.ext <- rstan::extract(object = glm)



  # get pmax
  glm.summary$pmax <- getPmax(glm.ext = glm.ext)




  # add gene id
  glm.summary$gene_name <- usage.data$gene_names



  # ppc
  ppc.data <- list(
    ppc.repertoire = getPpc(glm.ext = glm.ext,
                            usage.data = usage.data,
                            hdi.level = hdi.level),
    ppc.gene = getGroupStats(glm.ext = glm.ext,
                             usage.data = usage.data,
                             hdi.level = hdi.level))



  # frequentist tests, merge data
  t.test.stats <- getTTestStats(usage.data = usage.data)
  u.test.stats <- getManUStats(usage.data = usage.data)
  test.stats <- merge(x = t.test.stats, y = u.test.stats, by = "gene_name")


  # result
  result <- list(glm = glm,
                 glm.summary = glm.summary,
                 ppc.data = ppc.data,
                 usage.data = usage.data,
                 test.stats = test.stats)

  return (result)
}
