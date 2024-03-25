
# Description
# Posterior-predictive check for repertoires
get_ppc_rep <- function(glm, 
                        ud, 
                        hdi_lvl) {
  
  pars <- c("Yhat_rep", "Yhat_rep_prop")
  
  # summaries
  yhat <- summary(object = glm, pars = pars,
                  prob = c(0.5, (1-hdi_lvl)/2, 1-(1-hdi_lvl)/2))$summary
  yhat <- data.frame(yhat)
  colnames(yhat) <- c("mean", "se_mean", "sd", "median",
                      "L", "H", "Neff", "Rhat")
  yhat[, c("Rhat", "Neff")] <- NULL
  yhat$par <- rownames(yhat)
  yhat$par_name <- do.call(rbind, strsplit(x = yhat$par, split = "\\["))[, 1]
  par_i <- do.call(rbind, strsplit(x = yhat$par, split = "\\["))[, 2]
  par_i <- gsub(pattern = "\\]", replacement = '', x = par_i)
  par_i <- do.call(rbind, strsplit(x = par_i, split = ','))
  class(par_i) <- "numeric"
  yhat$G <- par_i[, 1]
  yhat$R <- par_i[, 2]
  rm(par_i)
  
  # split
  yhat_c <- yhat[yhat$par_name == pars[1], ] # count
  yhat_p <- yhat[yhat$par_name == pars[2], ] # prop
  rm(yhat)
  
  yhat_c$par <- NULL
  yhat_c$par_name <- NULL
  colnames(yhat_c)[seq_len(length.out = 6)] <- paste(
    "ppc", colnames(yhat_c)[seq_len(length.out = 6)],
    "count", sep = "_")
  
  yhat_p$par <- NULL
  yhat_p$par_name <- NULL
  colnames(yhat_p)[seq_len(length.out = 6)] <- paste(
    "ppc", colnames(yhat_p)[seq_len(length.out = 6)],
    "prop", sep = "_")
  
  yhat <- merge(x = yhat_c, y = yhat_p, by = c("G", "R"))
  rm(yhat_c, yhat_p)
  
  yhat$sample_id <- NA
  yhat$gene_name <- NA
  yhat$observed_count <- NA
  yhat$observed_prop <- NA
  
  for(i in seq_len(length.out = nrow(yhat))) {
    yhat$sample_id[i] <- ud$sample_names[yhat$R[i]]
    yhat$gene_name[i] <- ud$gene_names[yhat$G[i]]
    yhat$observed_count[i] <- ud$Y[yhat$G[i], yhat$R[i]]
    yhat$observed_prop[i] <- yhat$observed_count[i]/ud$N[yhat$R[i]]
  }
  
  if(ud$has_replicates) {
    t <- ud$proc_ud[, c("sample_id", "individual_id", "replicate_id", 
                        "condition", "gene_name")]
  } 
  else {
    t <- ud$proc_ud[, c("sample_id", "individual_id", 
                        "condition", "gene_name")]
  }
  return(merge(x = yhat, y = t, by = c("sample_id", "gene_name"), all.x = TRUE))
}



# Description
# Posterior-predictive for biological conditions
get_ppc_condition <- function(glm,
                              ud,
                              hdi_lvl) {
  
  # summaries
  yhat <- summary(object = glm, 
                  pars = c("Yhat_condition_prop"),
                  prob = c(0.5, (1-hdi_lvl)/2, 1-(1-hdi_lvl)/2))
  yhat <- yhat$summary
  yhat <- data.frame(yhat)
  colnames(yhat) <- c("ppc_mean_prop", "ppc_se_mean_prop", 
                      "ppc_sd_mean_prop", "ppc_median_prop", 
                      "ppc_L_prop", "ppc_H_prop", 
                      "Neff", "Rhat")
  yhat[, c("Rhat", "Neff")] <- NULL
  
  yhat$par <- rownames(yhat)
  yhat$par_name <- do.call(
    rbind, strsplit(x = yhat$par, split = "\\["))[, 1]
  par_i <- do.call(rbind, strsplit(x = yhat$par, split = "\\["))[,2]
  par_i <- gsub(pattern = "\\]", replacement = '', x = par_i)
  par_i <- do.call(rbind, strsplit(x = par_i, split = ','))
  class(par_i) <- "numeric"
  
  if(ud$N_condition==1) {
    yhat$condition_id <- 1
    yhat$gene_id <- par_i
  } else {
    yhat$condition_id <- par_i[, 1]
    yhat$gene_id <- par_i[, 2]
  }
  
  # condition map
  condition_map <- data.frame(condition_name = ud$condition_names,
                              condition_id = ud$condition_id)
  condition_map <- condition_map[duplicated(condition_map)==FALSE,]
  rownames(condition_map) <- condition_map$condition_id
  yhat$gene_name <- ud$gene_names[yhat$gene_id]
  yhat$condition <- condition_map[as.character(
    yhat$condition_id), "condition_name"]
  
  yhat$gene_id <- NULL
  yhat$condition_id <- NULL
  yhat$par <- NULL
  yhat$par_name <- NULL
  rownames(yhat) <- NULL
  return (yhat)
}
