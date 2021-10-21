
# Description:
# LOO = leave-one-out
# usage.data: 4 columns
#   * sample_id: char column
#   * condition: char column
#   * gene_name: char column
#   * gene_usage_count: num column
LOO <- function(usage.data,
                mcmc.warmup = 500,
                mcmc.steps = 1500,
                mcmc.chains = 4,
                mcmc.cores = 1,
                hdi.level = 0.95,
                adapt.delta = 0.95,
                max.treedepth = 12) {
  
  # before check convert summarized experiment object to data.frame
  if(inherits(x = usage.data, what = "SummarizedExperiment") == TRUE) {
    usage.data.raw <- usage.data
    usage.data <- convertSummarizedExperiment(usage.data.se = usage.data.raw)
  }
  
  # check inputs
  checkInput(usage.data = usage.data,
             mcmc.chains = as.integer(x = mcmc.chains),
             mcmc.cores = as.integer(x = mcmc.cores),
             mcmc.steps = as.integer(x = mcmc.steps),
             mcmc.warmup = as.integer(x = mcmc.warmup),
             hdi.level = hdi.level)
  
  # unique repertoire names
  Rs <- unique(usage.data$sample_id)
  
  # extra-stop condition
  if(length(Rs) <= 2) {
    stop("To perform LOO you need to provide as input at least 3 repertoires")
  }
  
  
  # compile model
  message("Compiling model ... \n")
  rstan::rstan_options(auto_write = TRUE)
  model.file <- system.file("extdata", "zibb.stan",
                            package = "IgGeneUsage")
  model <- rstan::stan_model(file = model.file, auto_write = TRUE)
  
  
  loo.out <- vector(mode = "list", length = length(Rs))
  names(loo.out) <- Rs
  for(r in 1:length(Rs)) {
    message(paste("LOO step: ", r, "\n", sep = ''))
    
    # here subset data
    temp.usage.data <- usage.data[usage.data$sample_id != Rs[r], ]
    
    # run DGU
    out <- LOO_DGU(usage.data = temp.usage.data,
                   mcmc.warmup = mcmc.warmup,
                   mcmc.steps = mcmc.steps,
                   mcmc.chains = mcmc.chains,
                   mcmc.cores = mcmc.cores,
                   hdi.level = hdi.level,
                   adapt.delta = adapt.delta,
                   max.treedepth = max.treedepth,
                   model = model)
    
    # collect results
    out$loo.id <- r
    loo.out[[Rs[r]]] <- out
  }
  
  # compact to data.frame and return
  loo.out <- do.call(rbind, loo.out)
  return (list(loo.summary = loo.out))
}



# Description:
# Light-weight version of DGU
LOO_DGU <- function(usage.data,
                mcmc.warmup = 500,
                mcmc.steps = 1500,
                mcmc.chains = 4,
                mcmc.cores = 1,
                hdi.level = 0.95,
                adapt.delta = 0.95,
                max.treedepth = 12,
                model) {
  
  
  # before check convert summarized experiment object to data.frame
  if(inherits(x = usage.data, what = "SummarizedExperiment") == TRUE) {
    usage.data.raw <- usage.data
    usage.data <- convertSummarizedExperiment(usage.data.se = usage.data.raw)
  }
  
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
  contrast <- paste(unique(usage.data$Xorg[usage.data$X == 1]),
                    " - ", unique(usage.data$Xorg[usage.data$X == -1]),
                    sep = '')
  
  
  # setup control list
  control.list <- list(adapt_delta = adapt.delta,
                       max_treedepth = max.treedepth)
  
  # stan sampling
  # monitor subset of parameters -> memory concern
  pars.relevant <- c("alpha_sigma", "beta_sigma",
                     "beta_gene_sigma", "phi",
                     "tau", "beta", "alpha_gene",
                     "beta_gene")
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
  message("Computing summaries ... \n")
  glm.summary <- rstan::summary(object = glm, digits = 4,
                                pars = "beta_gene",
                                prob = c(0.5, (1-hdi.level)/2,
                                         1-(1-hdi.level)/2))
  glm.summary <- glm.summary$summary
  glm.summary <- data.frame(glm.summary)
  colnames(glm.summary) <- c("es_mean", "es_mean_se",
                             "es_sd", "es_median",
                             "es_L", "es_H",
                             "Neff", "Rhat")
  glm.summary$contrast <- contrast
  
  
  # extract data
  message("Posterior extraction ... \n")
  glm.ext <- rstan::extract(object = glm)
  
  
  # get pmax
  message("Computing probability of DGU ... \n")
  glm.summary$pmax <- getPmax(glm.ext = glm.ext)
  
  
  # add gene id
  glm.summary$gene_name <- usage.data$gene_names
  
  
  return (glm.summary)
}


