
# Description:
# GU = gene usage
# ud: 4 columns
#   * sample_id: char column
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
  
  udr <- ud
  ud <- get_gu_usage(u = udr)
  
  # setup control list
  control_list <- list(adapt_delta = adapt_delta,
                       max_treedepth = max_treedepth)
  
  if(ud$N_group == 1) {
    pars <- get_pars(model = "GU")
    model <- stanmodels$gu
    glm <- rstan::sampling(object = model,
                           data = ud,
                           chains = mcmc_chains,
                           cores = mcmc_cores,
                           iter = mcmc_steps,
                           warmup = mcmc_warmup,
                           algorithm = "NUTS",
                           control = control_list,
                           pars = pars)
    
    message("Computing summaries ... \n")
    gu_summary <- get_gu_summary_univar(glm = glm, hdi_lvl = hdi_lvl, ud = ud)
    dgu_summary <- NA
    dgu_prob_summary <- NA
  } 
  else {
    pars <- get_pars(model = "DGU")
    model <- stanmodels$dgu
    glm <- rstan::sampling(object = model,
                           data = ud,
                           chains = mcmc_chains,
                           cores = mcmc_cores,
                           iter = mcmc_steps,
                           warmup = mcmc_warmup,
                           algorithm = "NUTS",
                           control = control_list,
                           pars = pars)
    
    message("Computing summaries ... \n")
    gu_summary <- get_gu_summary_anova(glm = glm, hdi_lvl = hdi_lvl, ud = ud)
    dgu_summary <- get_dgu_summary(glm = glm, hdi_lvl = hdi_lvl, ud = ud)
    dgu_prob_summary <- get_dgu_prob_summary(glm = glm, hdi_lvl = hdi_lvl, ud = ud)
  }
  
  # ppc
  message("Computing posterior predictions ... \n")
  ppc <- list(
    ppc_rep = get_ppc_rep(glm = glm, ud = ud, hdi_lvl = hdi_lvl),
    ppc_condition = get_ppc_condition(glm = glm, ud = ud, hdi_lvl = hdi_lvl))
  
  # result pack
  return (list(dgu_summary = dgu_summary,
               dgu_prob_summary = dgu_prob_summary,
               gu_summary = gu_summary,
               glm = glm,
               ppc = ppc,
               ud = ud))
}
