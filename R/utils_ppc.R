# Description
# Posterior-predictive check in repertoires, wrapper for un-/paired analysis
get_ppc_rep <- function(glm, ud, hdi_lvl, analysis_type) {
  if(analysis_type=="unpaired") {
    return(get_ppc_rep_u(glm = glm, ud = ud, hdi_lvl = hdi_lvl))
  }
  if(analysis_type=="paired") {
    return(get_ppc_rep_p(glm = glm, ud = ud, hdi_lvl = hdi_lvl))
  }
}

# Description
# Posterior-predictive check in repertoires, unpaired analysis
get_ppc_rep_u <- function(glm,
                          ud,
                          hdi_lvl) {
  
  pars <- c("Yhat", "Yhat_rep")
  
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
  colnames(yhat_c)[base::seq_len(length.out = 6)] <- paste(
    "ppc", colnames(yhat_c)[base::seq_len(length.out = 6)],
    "count", sep = "_")
  
  yhat_p$par <- NULL
  yhat_p$par_name <- NULL
  colnames(yhat_p)[base::seq_len(length.out = 6)] <- paste(
    "ppc", colnames(yhat_p)[base::seq_len(length.out = 6)],
    "prop", sep = "_")
  
  yhat <- merge(x = yhat_c, y = yhat_p, by = c("G", "R"))
  rm(yhat_c, yhat_p)
  
  yhat$sample_id <- NA
  yhat$gene_name <- NA
  yhat$observed_count <- NA
  yhat$observed_prop <- NA
  
  for(i in base::seq_len(length.out = nrow(yhat))) {
    yhat$sample_id[i] <- ud$sample_names[yhat$R[i]]
    yhat$gene_name[i] <- ud$gene_names[yhat$G[i]]
    yhat$observed_count[i] <- ud$Y[yhat$G[i], yhat$R[i]]
    yhat$observed_prop[i] <- yhat$observed_count[i]/ud$N[yhat$R[i]]
  }
  
  yhat <- base::merge(x = yhat, 
                      y = ud$proc_ud[, c("sample_id", "condition")], 
                      by = "sample_id", all.x = TRUE)
  return (yhat)
}

# Description
# Posterior-predictive check in repertoires, paired analysis
get_ppc_rep_p <- function(glm,
                          ud,
                          hdi_lvl) {
  
  get_ppc_single <- function(glm, 
                             ud, 
                             hdi_lvl, 
                             pars) {
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
    colnames(yhat_c)[base::seq_len(length.out = 6)] <- paste(
      "ppc", colnames(yhat_c)[base::seq_len(length.out = 6)],
      "count", sep = "_")
    
    yhat_p$par <- NULL
    yhat_p$par_name <- NULL
    colnames(yhat_p)[base::seq_len(length.out = 6)] <- paste(
      "ppc", colnames(yhat_p)[base::seq_len(length.out = 6)],
      "prop", sep = "_")
    
    yhat <- merge(x = yhat_c, y = yhat_p, by = c("G", "R"))
    rm(yhat_c, yhat_p)
    
    yhat$condition <- NA
    yhat$sample_name <- NA
    yhat$gene_name <- NA
    yhat$observed_count <- NA
    yhat$observed_prop <- NA
    
    
    for(i in base::seq_len(length.out = nrow(yhat))) {
      yhat$sample_name[i] <- ud$sample_names[yhat$R[i]]
      yhat$gene_name[i] <- ud$gene_names[yhat$G[i]]
      if(pars[1] == "Yhat_1") {
        yhat$observed_count[i] <- ud$Y_1[yhat$G[i], yhat$R[i]]
        yhat$observed_prop[i] <- yhat$observed_count[i]/ud$N_1[yhat$R[i]]
      } else {
        yhat$observed_count[i] <- ud$Y_2[yhat$G[i], yhat$R[i]]
        yhat$observed_prop[i] <- yhat$observed_count[i]/ud$N_2[yhat$R[i]]
      }
    }
    return (yhat)
  }
  
  y1 <- get_ppc_single(glm = glm, ud = ud, hdi_lvl = hdi_lvl, 
                       pars = c("Yhat_1", "Yhat_rep_1"))
  y1$group <- "Yhat_rep_1"
  y2 <- get_ppc_single(glm = glm, ud = ud, hdi_lvl = hdi_lvl, 
                       pars = c("Yhat_2", "Yhat_rep_2"))
  y2$group <- "Yhat_rep_2"
  return(rbind(y1, y2))
}

# Description
# Posterior-predictive in a biological condition
get_ppc_condition <- function(glm,
                              ud,
                              hdi_lvl,
                              analysis_type) {
  
  # summaries
  yhat <- summary(object = glm, 
                  pars = c("Yhat_condition"),
                  prob = c(0.5, (1-hdi_lvl)/2, 1-(1-hdi_lvl)/2))
  yhat <- yhat$summary
  yhat <- base::data.frame(yhat)
  base::colnames(yhat) <- c("ppc_mean_prop", "ppc_se_mean_prop", 
                            "ppc_sd_mean_prop", "ppc_median_prop", 
                            "ppc_L_prop", "ppc_H_prop", 
                            "Neff", "Rhat")
  yhat[, c("Rhat", "Neff")] <- NULL
  yhat$par <- base::rownames(yhat)
  yhat$par_name <- base::do.call(
    rbind, base::strsplit(x = yhat$par, split = "\\["))[, 1]
  par_i <- do.call(rbind, base::strsplit(x = yhat$par, split = "\\["))[, 2]
  par_i <- base::gsub(pattern = "\\]", replacement = '', x = par_i)
  par_i <- base::do.call(rbind, base::strsplit(x = par_i, split = ','))
  base::class(par_i) <- "numeric"
  yhat$X <- par_i[, 1]
  yhat$G <- par_i[, 2]
  yhat$gene_name <- ud$gene_names[yhat$G]
  if(analysis_type=="unpaired") {
    yhat$condition <- ifelse(test = yhat$X == 1, 
                             yes = ud$X_org[ud$X==1][1],
                             no = ud$X_org[ud$X==-1][1])
  }
  
  return (yhat)
}