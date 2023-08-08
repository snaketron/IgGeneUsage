
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
  
  # check inputs
  check_dgu_input(ud = ud,
                  mcmc_chains = base::as.integer(x = mcmc_chains),
                  mcmc_cores = base::as.integer(x = mcmc_cores),
                  mcmc_steps = base::as.integer(x = mcmc_steps),
                  mcmc_warmup = base::as.integer(x = mcmc_warmup),
                  hdi_lvl = hdi_lvl)
  
  # browser()
  # format input usage
  udr <- ud
  ud <- get_usage(u = udr)
  
  analysis <- get_analysis_type(udr)
  
  # setup control list
  control_list <- base::list(adapt_delta = adapt_delta,
                             max_treedepth = max_treedepth)
  
  # stan sampling
  if(analysis == "paired") {
    # model <- rstan::stan_model(file = "inst/stan/dgu_pair.stan")
    model <- stanmodels$dgu_pair
    pars <- get_pars(model = "DGU", analysis = "paired")
  } 
  else {
    # model <- rstan::stan_model(file = "inst/stan/dgu_unpair.stan")
    model <-  stanmodels$dgu_unpair
    pars <- get_pars(model = "DGU", analysis = "unpaired")
  }
  glm <- rstan::sampling(object = model,
                         data = ud,
                         chains = mcmc_chains,
                         cores = mcmc_cores,
                         iter = mcmc_steps,
                         warmup = mcmc_warmup,
                         algorithm = "NUTS",
                         control = control_list,
                         pars = pars)
  
  
  # get summary
  message("Computing summaries ... \n")
  glm_summary <- get_glm_summary_dgu(glm = glm, 
                                     hdi_lvl = hdi_lvl, 
                                     ud = ud)
  

  # ppc
  message("Computing posterior predictions ... \n")
  ppc <- list(
    ppc_rep = get_ppc_rep(glm = glm, 
                          ud = ud, 
                          hdi_lvl = hdi_lvl, 
                          analysis = analysis),
    ppc_condition = get_ppc_condition(glm = glm, 
                                      ud = ud, 
                                      hdi_lvl = hdi_lvl,
                                      analysis = analysis))
  
  # frequentist tests, merge data
  message("Computing frequentist DGU ... \n")
  t_test_stats <- get_ttest(ud = ud, paired = analysis == "paired")
  u_test_stats <- get_manu(ud = ud, paired = analysis == "paired")
  test_summary <- merge(x = t_test_stats, y = u_test_stats, by = "gene_name")
  
  # result pack
  return (list(glm_summary = glm_summary,
               test_summary = test_summary,
               glm = glm,
               ppc = ppc,
               ud = ud))
}
