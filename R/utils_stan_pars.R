# Which model parameters should be monitored?
# model = DGU or GU_univar or GU_anova
get_pars <- function(model) {
  if(model == "GU") {
    return(c("phi",
             "kappa",
             "alpha",
             "alpha_gene_mu",
             "alpha_gene_sigma",
             # "beta_gene_mu",
             #"beta_pop_sigma", 
             # "beta_gene_sigma",
             "log_lik",
             "Yhat",
             "Yhat_rep",
             "Yhat_condition",
             "theta_condition",
             "theta_rep"))
             #"beta_alpha",
             #"beta_nu",
             #"beta_tau"))
  }
  
  if(model == "DGU") {
    return(c("beta",
             "beta_pop_sigma", 
             "beta_gene_sigma",
             "phi",
             "kappa",
             "alpha_gene_mu", 
             "beta_gene_mu",
             "log_lik",
             "Yhat",
             "Yhat_rep",
             "Yhat_condition",
             "dgu",
             "dgu_prob",
             "theta_condition",
             "theta_rep"))
  }
  stop("wrong stan pars")
}