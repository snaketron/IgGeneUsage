
get_gu_summary_anova <- function(glm, hdi_lvl, ud) {
  
  gu_summary <- rstan::summary(object = glm, 
                               digits = 4,
                               pars = "theta_condition",
                               prob = c(0.5, (1-hdi_lvl)/2,
                                        1-(1-hdi_lvl)/2))
  
  gu_summary <- gu_summary$summary
  gu_summary <- data.frame(gu_summary)
  colnames(gu_summary) <- c("prob_mean", "prob_mean_se",
                            "prob_sd", "prob_median",
                            "prob_L", "prob_H",
                            "Neff", "Rhat")
  gu_summary[, c("Rhat", "Neff")] <- NULL
  
  par <- rownames(gu_summary)
  par <- gsub(pattern = "theta_condition|\\[|\\]", replacement = '', x = par)
  par <- do.call(rbind, strsplit(x = par, split = ','))
  
  gu_summary$gene_id <- as.numeric(par[,2])
  gu_summary$gene_name <- ud$gene_names[gu_summary$gene_id]
  gu_summary$condition_id <- as.numeric(par[,1])
  gu_summary$condition <- NA
  for(i in 1:ud$N_group) {
    gu_summary$condition[which(gu_summary$condition_id == i)]<-
      ud$group_names[which(ud$group_id == i)[1]]
  }
  
  # remove unused vars
  gu_summary$gene_id <- NULL
  gu_summary$condition_id <- NULL
  rownames(gu_summary) <- NULL
  return(gu_summary)
}


get_gu_summary_univar <- function(glm, hdi_lvl, ud) {
  
  gu_summary <- rstan::summary(object = glm, 
                               digits = 4,
                               pars = "theta_condition",
                               prob = c(0.5, (1-hdi_lvl)/2,
                                        1-(1-hdi_lvl)/2))
  
  gu_summary <- gu_summary$summary
  gu_summary <- data.frame(gu_summary)
  colnames(gu_summary) <- c("prob_mean", "prob_mean_se",
                            "prob_sd", "prob_median",
                            "prob_L", "prob_H",
                            "Neff", "Rhat")
  gu_summary[, c("Rhat", "Neff")] <- NULL
  
  par <- rownames(gu_summary)
  par <- gsub(pattern = "theta_condition|\\[|\\]", replacement = '', x = par)
  par <- do.call(rbind, strsplit(x = par, split = ','))
  
  gu_summary$gene_id <- as.numeric(par)
  gu_summary$gene_name <- ud$gene_names[gu_summary$gene_id]
  gu_summary$condition <- ud$group_names[1]
  
  # remove unused vars
  gu_summary$gene_id <- NULL
  rownames(gu_summary) <- NULL
  return(gu_summary)
}

