
# Description:
# DGU = differential gene usage
# usage.data: 4 columns
#   * sample_id: char column
#   * condition: char column
#   * gene_name: char column
#   * gene_usage_count: num column
DGU <- function(usage.data,
                mcmc.warmup = 500,
                mcmc.steps = 1500,
                mcmc.chains = 4,
                mcmc.cores = 4,
                hdi.level = 0.95,
                adapt.delta = 0.95,
                max.treedepth = 12) {


  # check inputs
  checkInput(usage.data = usage.data,
             mcmc.chains = as.integer(x = mcmc.chains),
             mcmc.cores = as.integer(x = mcmc.cores),
             mcmc.steps = as.integer(x = mcmc.steps),
             mcmc.warmup = as.integer(x = mcmc.warmup),
             hdi.level = hdi.level)


  # format input usage
  usage.data.raw <- usage.data
  usage.data <- getUsageData(usage = usage.data.raw)


  # contrast
  contrast <- paste("Contrast: ", unique(usage.data$Xorg[usage.data$X == 1]),
                    " - ", unique(usage.data$Xorg[usage.data$X == -1]),
                    sep = '')


  # model
  cat("Compiling model ... \n")
  rstan::rstan_options(auto_write = TRUE)
  model.file <- system.file("extdata", "zibb.stan",
                            package = "IgGeneUsage")
  model <- rstan::stan_model(file = model.file,
                             auto_write = TRUE)


  # setup control list
  control.list <- list(adapt_delta = adapt.delta,
                       max_treedepth = max.treedepth)


  # stan sampling
  # monitor subset of parameters -> memory concern
  pars.relevant <- c("alpha_grand", "beta_grand",
                     "alpha_sigma", "beta_sigma",
                     "beta_gene_sigma", "phi",
                     "tau", "beta", "alpha_gene",
                     "beta_gene", "log_lik",
                     "Yhat",  "Yhat_individual",
                     "Yhat_gene")
  glm <- rstan::sampling(object = model,
                         data = usage.data,
                         chains = mcmc.chains,
                         cores = mcmc.cores,
                         iter = mcmc.steps,
                         warmup = mcmc.warmup,
                         refresh = 250,
                         control = control.list,
                         pars = pars.relevant,
                         algorithm = "NUTS")


  # get summary
  cat("Computing summaries ... \n")
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
  cat("Posterior extraction ... \n")
  glm.ext <- rstan::extract(object = glm)


  # get pmax
  cat("Computing probability of DGU ... \n")
  glm.summary$pmax <- getPmax(glm.ext = glm.ext)


  # add gene id
  glm.summary$gene_name <- usage.data$gene_names


  # ppc
  cat("Computing posterior predictions ... \n")
  ppc.data <- list(
    ppc.repertoire = getPpc(glm.ext = glm.ext,
                            usage.data = usage.data,
                            hdi.level = hdi.level),
    ppc.gene = getGroupStats(glm.ext = glm.ext,
                             usage.data = usage.data,
                             hdi.level = hdi.level))



  # frequentist tests, merge data
  cat("Computing frequentist DGU ... \n")
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
