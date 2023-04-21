
# Description:
# DGU = differential gene usage
# ud: 4 columns
#   * sample_id: char column
#   * condition: char column
#   * gene_name: char column
#   * gene_usage_count: num column
DGU <- function(ud,
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
    ud <- get_SummarizedExperiment(ud_se = udr)
  }
  
  # check inputs
  check_input(ud = ud,
              mcmc_chains = base::as.integer(x = mcmc_chains),
              mcmc_cores = base::as.integer(x = mcmc_cores),
              mcmc_steps = base::as.integer(x = mcmc_steps),
              mcmc_warmup = base::as.integer(x = mcmc_warmup),
              hdi_lvl = hdi_lvl)
  
  # format input usage
  udr <- ud
  ud <- get_usage(u = udr)
  
  # contrast
  contrast <- paste0(unique(ud$Xorg[ud$X == 1]), " - ", 
                     unique(ud$Xorg[ud$X == -1]))
  
  # setup control list
  control_list <- list(adapt_delta = adapt_delta,
                       max_treedepth = max_treedepth)
  
  # stan sampling
  glm <- rstan::sampling(object = stanmodels$zibb,
                         data = ud,
                         chains = mcmc_chains,
                         cores = mcmc_cores,
                         iter = mcmc_steps,
                         warmup = mcmc_warmup,
                         algorithm = "NUTS",
                         control = control_list,
                         pars = c("alpha_sigma", "beta_sigma",
                                  "beta_gene_sigma", "phi",
                                  "tau", 
                                  "z",
                                  "alpha_gene", "beta_gene", 
                                  "log_lik", "Yhat",  
                                  "Yhat_individual", "Yhat_gene"))
  
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
  glm_summary[, c("Rhat", "Neff")] <- NULL
  glm_summary$contrast <- contrast
  
  # extract data
  message("Posterior extraction ... \n")
  glm_ext <- rstan::extract(object = glm)
  
  # get pmax
  message("Computing probability of DGU ... \n")
  glm_summary$pmax <- get_pmax(glm_ext = glm_ext)
  
  # add gene id
  glm_summary$gene_name <- ud$gene_names
  
  # ppc
  message("Computing posterior predictions ... \n")
  ppc <- list(
    ppc_repertoire = get_ppc_rep(glm = glm, ud = ud, hdi_lvl = hdi_lvl),
    ppc_gene = get_ppc_gene(glm = glm, ud = ud, hdi_lvl = hdi_lvl))
  
  # frequentist tests, merge data
  message("Computing frequentist DGU ... \n")
  t_test_stats <- get_ttest(ud = ud)
  u_test_stats <- get_manu(ud = ud)
  test_summary <- merge(x = t_test_stats, y = u_test_stats, by = "gene_name")
  
  # result pack
  return (list(glm_summary = glm_summary,
               test_summary = test_summary,
               glm = glm,
               ppc_data = ppc,
               ud = ud))
}
