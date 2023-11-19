# Which model parameters should be monitored?
# model = DGU or GU_univar or GU_anova
get_pars <- function(model) {
  if(model == "GU") {
    return(c("phi",
             "kappa",
             "alpha_gene_mu",
             "alpha_gene_sigma",
             "theta",
             "Yhat_rep",
             "Yhat_rep_prop",
             "Yhat_condition_prop",
             "log_lik"))
  }
  
  if(model == "DGU") {
    return(c("phi",
             "kappa",
             "alpha_gene_mu",
             "beta_gene_mu",
             "beta_gene_sigma",
             "beta_pop_sigma",
             "theta",
             "Yhat_rep",
             "Yhat_rep_prop",
             "Yhat_condition_prop",
             "log_lik",
             "dgu",
             "dgu_prob"))
  }
  stop("wrong stan pars")
}