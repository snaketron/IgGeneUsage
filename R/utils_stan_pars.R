# Which model parameters should be monitored?
# model = DGU or GU_univar or GU_anova
# analysis = unpaired or paired (only relevant if model = DGU)
get_pars <- function(model, analysis) {
  if(model == "GU") {
    return(c("alpha_pop_mu", 
             "alpha_pop_sigma",
             "phi",
             "z",
             "alpha_gene_mu",
             "log_lik",
             "Yhat",
             "Yhat_rep",
             "Yhat_condition",
             "prob_gene"))
  }
  
  if(model == "DGU") {
    return(c("beta",
             "alpha_pop_mu",
             "alpha_pop_sigma", 
             "beta_pop_sigma", 
             "beta_gene_sigma",
             "phi",
             "z",
             "alpha_gene_mu", 
             "beta_gene_mu",
             "log_lik",
             "Yhat",
             "Yhat_rep",
             "Yhat_condition",
             "dgu"))
  }
  stop("wrong stan pars")
}