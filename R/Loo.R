
# Description:
# LOO = leave-one-out
# ud: 4 columns
#   * sample_id: char column
#   * condition: char column
#   * gene_name: char column
#   * gene_usage_count: num column
LOO <- function(ud,
                mcmc_warmup = 500,
                mcmc_steps = 1500,
                mcmc_chains = 4,
                mcmc_cores = 1,
                hdi_lvl = 0.95,
                adapt_delta = 0.95,
                max_treedepth = 12) {
  
  # before check convert summarized experiment object to data.frame
  if(inherits(x = ud, what = "SummarizedExperiment") == TRUE) {
    udr <- ud
    ud <- convertSummarizedExperiment(ud_se = udr)
  }
  
  # check inputs
  checkInput(ud = ud,
             mcmc_chains = as.integer(x = mcmc_chains),
             mcmc_cores = as.integer(x = mcmc_cores),
             mcmc_steps = as.integer(x = mcmc_steps),
             mcmc_warmup = as.integer(x = mcmc_warmup),
             hdi_lvl = hdi_lvl)
  
  # unique repertoire names
  Rs <- unique(ud$sample_id)
  
  # extra-stop condition
  if(length(Rs) <= 2) {
    stop("To perform LOO you need to provide as input at least 3 repertoires")
  }
  
  loo_out <- vector(mode = "list", length = length(Rs))
  names(loo_out) <- Rs
  for(r in base::seq_len(length.out = length(Rs))) {
    base::message("LOO step: ", r, "\n", sep = '')
    
    # here subset data
    temp_ud <- ud[ud$sample_id != Rs[r], ]
    
    # run DGU
    out <- LOO_DGU(ud = temp_ud,
                   mcmc_warmup = mcmc_warmup,
                   mcmc_steps = mcmc_steps,
                   mcmc_chains = mcmc_chains,
                   mcmc_cores = mcmc_cores,
                   hdi_lvl = hdi_lvl,
                   adapt_delta = adapt_delta,
                   max_treedepth = max_treedepth)
    
    # collect results
    out$loo.id <- r
    loo_out[[Rs[r]]] <- out
  }
  
  # compact to data.frame and return
  loo_out <- do.call(rbind, loo_out)
  return (list(loo.summary = loo_out))
}



# Description:
# Light-weight version of DGU
LOO_DGU <- function(ud,
                    mcmc_warmup = 500,
                    mcmc_steps = 1500,
                    mcmc_chains = 4,
                    mcmc_cores = 1,
                    hdi_lvl = 0.95,
                    adapt_delta = 0.95,
                    max_treedepth = 12) {
  
  
  # before check convert summarized experiment object to data.frame
  if(inherits(x = ud, what = "SummarizedExperiment") == TRUE) {
    udr <- ud
    ud <- convertSummarizedExperiment(ud_se = udr)
  }
  
  # check inputs
  checkInput(ud = ud,
             mcmc_chains = base::as.integer(x = mcmc_chains),
             mcmc_cores = base::as.integer(x = mcmc_cores),
             mcmc_steps = base::as.integer(x = mcmc_steps),
             mcmc_warmup = base::as.integer(x = mcmc_warmup),
             hdi_lvl = hdi_lvl)
  
  # format input usage
  udr <- ud
  ud <- get_usage(u = udr)
  
  # contrast
  contrast <- paste0(unique(ud$Xorg[ud$X == 1]),
                     " - ", 
                     unique(ud$Xorg[ud$X == -1]))
  
  # setup control list
  control_list <- list(adapt_delta = adapt_delta,
                       max_treedepth = max_treedepth)
  
  # stan sampling
  # monitor subset of parameters -> memory concern
  pars_rel <- c("alpha_sigma", "beta_sigma",
                "beta_gene_sigma", "phi",
                "tau", "beta", "alpha_gene",
                "beta_gene")
  glm <- rstan::sampling(object = stanmodels$zibb,
                         data = ud,
                         chains = mcmc_chains,
                         cores = mcmc_cores,
                         iter = mcmc_steps,
                         warmup = mcmc_warmup,
                         refresh = 250,
                         control = control_list,
                         pars = pars_rel,
                         algorithm = "NUTS")
  
  # get summary
  message("Computing summaries ... \n")
  glm_summary <- rstan::summary(object = glm, 
                                digits = 4,
                                pars = "beta_gene",
                                prob = c(0.5, (1-hdi_lvl)/2,
                                         1-(1-hdi_lvl)/2))
  glm_summary <- glm_summary$summary
  glm_summary <- data.frame(glm_summary)
  colnames(glm_summary) <- c("es_mean", "es_mean_se",
                             "es_sd", "es_median",
                             "es_L", "es_H",
                             "Neff", "Rhat")
  glm_summary$contrast <- contrast
  
  # extract data
  message("Posterior extraction ... \n")
  glm_ext <- rstan::extract(object = glm)
  
  # get pmax
  message("Computing probability of DGU ... \n")
  glm_summary$pmax <- getPmax(glm_ext = glm_ext)
  
  # add gene id
  glm_summary$gene_name <- ud$gene_names
  
  return(glm_summary)
}


