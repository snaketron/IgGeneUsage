
get_condition_prop <- function(glm, hdi_lvl, ud, model_type) {
  if(model_type == "DGU") {
    return(get_condition_prop_dgu(glm = glm, hdi_lvl = hdi_lvl, ud = ud))
  }
  if(model_type == "GU") {
    return(get_condition_prop_gu(glm = glm, hdi_lvl = hdi_lvl, ud = ud))
  }
}


get_condition_prop_dgu <- function(glm, hdi_lvl, ud) {
  
  gu_summary <- summary(object = glm, 
                        digits = 4,
                        pars = "Yhat_condition_prop",
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
  
  par <- gsub(pattern = "Yhat_condition_prop|\\[|\\]", 
              replacement = '', x = par)
  par <- do.call(rbind, strsplit(x = par, split = ','))
  
  gu_summary$gene_id <- as.numeric(par[,2])
  gu_summary$gene_name <- ud$gene_names[gu_summary$gene_id]
  gu_summary$condition_id <- as.numeric(par[,1])
  gu_summary$condition <- NA
  for(i in 1:ud$N_condition) {
    gu_summary$condition[which(gu_summary$condition_id == i)]<-
      ud$condition_names[which(ud$condition_id == i)[1]]
  }
  
  # remove unused vars
  gu_summary$gene_id <- NULL
  gu_summary$condition_id <- NULL
  rownames(gu_summary) <- NULL
  return(gu_summary)
}


get_condition_prop_gu <- function(glm, hdi_lvl, ud) {
  
  gu_summary <- summary(object = glm, 
                        digits = 4,
                        pars = "Yhat_condition_prop",
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
  par <- gsub(pattern = "Yhat_condition_prop|\\[|\\]", 
              replacement = '', x = par)
  par <- do.call(rbind, strsplit(x = par, split = ','))
  
  gu_summary$gene_id <- as.numeric(par)
  gu_summary$gene_name <- ud$gene_names[gu_summary$gene_id]
  gu_summary$condition <- ud$condition_names[1]
  
  # remove unused vars
  gu_summary$gene_id <- NULL
  rownames(gu_summary) <- NULL
  return(gu_summary)
}


get_sample_prop_gu <- function(glm, hdi_lvl, ud) {
  
  gu <- summary(object = glm, digits = 4, pars = "theta",
                prob = c(0.5, (1-hdi_lvl)/2, 1-(1-hdi_lvl)/2))
  gu <- data.frame(gu$summary)
  colnames(gu) <- c("theta_mean", "theta_mean_se",
                    "theta_sd", "theta_median",
                    "theta_L", "theta_H",
                    "Neff", "Rhat")
  gu[, c("Rhat", "Neff")] <- NULL
  
  par <- rownames(gu)
  par <- gsub(pattern = "theta|\\[|\\]", replacement = '', x = par)
  par <- do.call(rbind, strsplit(x = par, split = ','))
  
  gu$gene_id <- as.numeric(par[,2])
  gu$sample_id <- as.numeric(par[,1])
  
  gu$gene_name <- ud$gene_names[gu$gene_id]
  gu$sample_name <- ud$sample_names[gu$sample_id]
  
  m <- ud$proc_ud[, c("sample_id", "individual_id", "individual_org_name",
                      "gene_name", "condition", "gene_usage_prop")]
  
  gu <- merge(x = gu, y = m, 
              by.x = c("sample_id", "gene_name"), 
              by.y = c("sample_id", "gene_name"))
  
  # remove unused vars
  gu$gene_id <- NULL
  gu$sample_id <- NULL
  rownames(gu) <- NULL
  return(gu)
}

