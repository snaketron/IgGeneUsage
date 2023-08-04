
# Description:
# GU = gene usage
# ud: 4 columns
#   * sample_id: char column
#   * gene_name: char column
#   * gene_usage_count: num column
GU <- function(ud,
               mcmc_warmup = 500,
               mcmc_steps = 1500,
               mcmc_chains = 4,
               mcmc_cores = 1,
               hdi_lvl = 0.95,
               adapt_delta = 0.95,
               max_treedepth = 12) {
  
  # check inputs
  # check_gu_input(ud = ud,
  #                mcmc_chains = base::as.integer(x = mcmc_chains),
  #                mcmc_cores = base::as.integer(x = mcmc_cores),
  #                mcmc_steps = base::as.integer(x = mcmc_steps),
  #                mcmc_warmup = base::as.integer(x = mcmc_warmup),
  #                hdi_lvl = hdi_lvl)
  analysis_type <- "unpaired"
  udr <- ud
  ud <- get_gu_usage(u = udr)
  
  # setup control list
  control_list <- list(adapt_delta = adapt_delta,
                       max_treedepth = max_treedepth)
  
  
  if(ud$N_group == 1) {
    
    pars <- c("alpha_pop_mu",
              "alpha_pop_sigma",
              "alpha_gene_sigma",
              "phi",
              "z", "z_mu", "z_phi",
              "alpha_gene_mu",
              "log_lik",
              "Yhat",
              "Yhat_rep",
              "Yhat_condition",
              "prob_gene")
    
    m_gu_univar <- rstan::stan_model(file = "inst/stan/gu_univar.stan")
    glm <- rstan::sampling(object = m_gu_univar,#stanmodels$gu_univar,
                           data = ud,
                           chains = mcmc_chains,
                           cores = mcmc_cores,
                           iter = mcmc_steps,
                           warmup = mcmc_warmup,
                           algorithm = "NUTS",
                           control = control_list,
                           pars = pars)
    
    message("Computing summaries ... \n")
    glm_summary <- get_glm_summary_gu_univar(glm = glm, 
                                             hdi_lvl = hdi_lvl, 
                                             ud = ud)
  } 
  else {
    pars <- c("beta",
              "alpha_pop_mu",
              "alpha_pop_sigma", "beta_pop_sigma",
              "alpha_gene_sigma", "beta_gene_sigma",
              "phi",
              "z", "z_mu", "z_phi",
              "alpha_gene_mu", "beta_gene_mu",
              "log_lik",
              "Yhat",
              "Yhat_rep",
              "Yhat_condition")
    
    m_gu_anova <- rstan::stan_model(file = "inst/stan/gu_anova.stan")
    glm <- rstan::sampling(object = m_gu_anova,#stanmodels$gu_anova,
                           data = ud,
                           chains = mcmc_chains,
                           cores = mcmc_cores,
                           iter = mcmc_steps,
                           warmup = mcmc_warmup,
                           algorithm = "NUTS",
                           control = control_list,
                           pars = pars)
    
    message("Computing summaries ... \n")
    glm_summary <- get_glm_summary_gu_anova(glm = glm, 
                                            hdi_lvl = hdi_lvl, 
                                            ud = ud)
  }
  
  # ppc
  message("Computing posterior predictions ... \n")
  ppc <- list(
    ppc_rep = get_ppc_rep(glm = glm,
                          ud = ud,
                          hdi_lvl = hdi_lvl,
                          analysis_type = analysis_type),
    ppc_condition = get_ppc_condition(glm = glm,
                                      ud = ud,
                                      hdi_lvl = hdi_lvl,
                                      analysis_type = analysis_type))
  
  # result pack
  return (list(glm_summary = glm_summary,
               glm = glm,
               ud = ud))
}
