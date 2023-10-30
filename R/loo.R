
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
  
  # check inputs
  check_dgu_input(ud = ud,
                  mcmc_chains = base::as.integer(x = mcmc_chains),
                  mcmc_cores = base::as.integer(x = mcmc_cores),
                  mcmc_steps = base::as.integer(x = mcmc_steps),
                  mcmc_warmup = base::as.integer(x = mcmc_warmup),
                  hdi_lvl = hdi_lvl)
  
  # process data
  ud <- get_gu_usage(u = ud)
  ud <- ud$proc_ud[, c("sample_id","condition","gene_name","gene_usage_count")]
  
  # setup control list
  control_list <- list(adapt_delta = adapt_delta,
                       max_treedepth = max_treedepth)
  
  # unique repertoire names
  ud$loo_id <- ud$sample_id
  rs <- unique(ud$loo_id)
  
  # extra-stop condition
  if(length(rs) <= 2) {
    stop("To perform LOO you need to provide as input at least 3 repertoires")
  }
  
  loo_out <- vector(mode = "list", length = length(rs))
  names(loo_out) <- rs
  for(r in base::seq_len(length.out = length(rs))) {
    base::message("LOO step: ", r, "\n", sep = '')
    
    # here subset data
    temp_ud <- ud[ud$loo_id != rs[r], ]
    temp_ud$loo_id <- NULL
    
    # run DGU
    out <- DGU(ud = temp_ud,
               mcmc_warmup = mcmc_warmup,
               mcmc_steps = mcmc_steps,
               mcmc_chains = mcmc_chains,
               mcmc_cores = mcmc_cores,
               hdi_lvl = hdi_lvl,
               adapt_delta = adapt_delta,
               max_treedepth = max_treedepth)
    
    if(is.data.frame(out$gu_summary)==TRUE) {
      out$gu_summary$loo_id <- rs[r]
    }
    if(is.data.frame(out$dgu_summary)==TRUE) {
      out$dgu_summary$loo_id <- rs[r]
    }
    if(is.data.frame(out$dgu_prob_summary)==TRUE) {
      out$dgu_prob_summary$loo_id <- rs[r]
    }
    
    # collect results
    loo_out[[rs[r]]] <- out
  }
  
  return (loo_out)
}


