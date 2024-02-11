
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
  ud <- get_usage(u = udr)
  
  # setup control list
  control_list <- list(adapt_delta = adapt_delta,
                       max_treedepth = max_treedepth)
  
  # get model
  m <- get_model(has_conditions = ud$has_conditions, 
                 has_replicates = ud$has_replicates)
  
  # fit model
  glm <- rstan::sampling(object = m$model,
                         data = ud,
                         chains = mcmc_chains,
                         cores = mcmc_cores,
                         iter = mcmc_steps,
                         warmup = mcmc_warmup,
                         algorithm = "NUTS",
                         control = control_list,
                         pars = m$pars,
                         refresh = 50)
  
  if(m$model_type=="GU") {
    message("Computing summaries ... \n")
    gu <- get_gu_summary_gu(glm = glm, hdi_lvl = hdi_lvl, ud = ud)
    dgu <- NA
    dgu_prob <- NA
  }
  if(m$model_type=="DGU") {
    message("Computing summaries ... \n")
    gu <- get_gu_summary_dgu(glm = glm, hdi_lvl = hdi_lvl, ud = ud)
    dgu <- get_dgu_summary(glm = glm, hdi_lvl = hdi_lvl, ud = ud)
    dgu_prob <- get_dgu_prob_summary(glm = glm, hdi_lvl = hdi_lvl, ud = ud)
  }
  
  # ppc
  message("Computing posterior predictions ... \n")
  ppc <- list(
    ppc_rep = get_ppc_rep(glm = glm, ud = ud, hdi_lvl = hdi_lvl),
    ppc_condition = get_ppc_condition(glm = glm, ud = ud, hdi_lvl = hdi_lvl))
  
  # result pack
  return (list(dgu = dgu,
               dgu_prob = dgu_prob,
               gu = gu,
               glm = glm,
               ppc = ppc,
               ud = ud))
}
