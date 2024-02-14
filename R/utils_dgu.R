get_contrast_map <- function(ud) {
  # condition map
  condition_map <- data.frame(condition_name = ud$condition_names,
                              condition_id = ud$condition_id)
  condition_map <- condition_map[duplicated(condition_map)==FALSE,]
  rownames(condition_map) <- condition_map$condition_id
  
  # contrast map
  c_n <- max(condition_map$condition_id)*(max(condition_map$condition_id)-1)/2
  c_map <- data.frame(c_id = seq_len(c_n))
  c_names <- c()
  for(i in 1:(max(condition_map$condition_id)-1)) {
    for(j in (i+1):max(condition_map$condition_id)) {
      c_names <- c(c_names,
                   paste0(condition_map[as.character(i), "condition_name"],
                          "-vs-",
                          condition_map[as.character(j), "condition_name"]))
    }
  }
  c_map$contrast <- c_names
  return(c_map)
}


get_dgu <- function(glm, hdi_lvl, ud, model_type) {
  if(model_type == "GU") {
    return(NA)
  }
  contrast_map <- get_contrast_map(ud)
  
  dgu_summary <- summary(object = glm, 
                         digits = 4,
                         pars = "dgu",
                         prob = c(0.5, (1-hdi_lvl)/2,
                                  1-(1-hdi_lvl)/2))
  
  dgu_summary <- dgu_summary$summary
  dgu_summary <- data.frame(dgu_summary)
  colnames(dgu_summary) <- c("es_mean", "es_mean_se",
                             "es_sd", "es_median",
                             "es_L", "es_H",
                             "Neff", "Rhat")
  dgu_summary[, c("Rhat", "Neff")] <- NULL
  
  par <- rownames(dgu_summary)
  par <- gsub(pattern = "dgu|\\[|\\]", replacement = '', x = par)
  par <- do.call(rbind, strsplit(x = par, split = ','))
  
  
  dgu_summary$gene_id <- as.numeric(par[,1])
  dgu_summary$c_id <- as.numeric(par[,2])
  dgu_summary <- merge(x = dgu_summary, 
                       y = contrast_map, 
                       by = "c_id", 
                       all.x = TRUE)
  
  dgu_summary$gene_name <- ud$gene_names[dgu_summary$gene_id]
  dgu_summary$pmax <- NA
  
  # get pmax
  glm_ext <- rstan::extract(object = glm, par = "dgu")$dgu
  for(i in seq_len(nrow(dgu_summary))) {
    dgu_summary$pmax[i] <- get_pmax(
      x = glm_ext[,dgu_summary$gene_id[i], dgu_summary$c_id[i]])
  }
  
  # remove unused vars
  dgu_summary$gene_id <- NULL
  dgu_summary$c_id <- NULL
  rownames(dgu_summary) <- NULL
  return(dgu_summary)
}


get_dgu_prob <- function(glm, hdi_lvl, ud, model_type) {
  if(model_type == "GU") {
    return(NA)
  }
  contrast_map <- get_contrast_map(ud)
  
  dgu_summary <- summary(object = glm, 
                         digits = 4,
                         pars = "dgu_prob",
                         prob = c(0.5, (1-hdi_lvl)/2,
                                  1-(1-hdi_lvl)/2))
  
  dgu_summary <- dgu_summary$summary
  dgu_summary <- data.frame(dgu_summary)
  colnames(dgu_summary) <- c("es_prob_mean", "es_prob_mean_se",
                             "es_prob_sd", "es_prob_median",
                             "es_prob_L", "es_prob_H",
                             "Neff", "Rhat")
  dgu_summary[, c("Rhat", "Neff")] <- NULL
  
  par <- rownames(dgu_summary)
  par <- gsub(pattern = "dgu_prob|\\[|\\]", replacement = '', x = par)
  par <- do.call(rbind, strsplit(x = par, split = ','))
  
  
  dgu_summary$gene_id <- as.numeric(par[,1])
  dgu_summary$c_id <- as.numeric(par[,2])
  dgu_summary <- merge(x = dgu_summary, 
                       y = contrast_map, 
                       by = "c_id", 
                       all.x = TRUE)
  
  dgu_summary$gene_name <- ud$gene_names[dgu_summary$gene_id]
  dgu_summary$pmax <- NA
  
  # get pmax
  glm_ext <- rstan::extract(object = glm, par = "dgu_prob")$dgu_prob
  for(i in seq_len(nrow(dgu_summary))) {
    dgu_summary$pmax[i] <- get_pmax(
      x = glm_ext[,dgu_summary$gene_id[i], dgu_summary$c_id[i]])
  }
  
  # remove unused vars
  dgu_summary$gene_id <- NULL
  dgu_summary$c_id <- NULL
  rownames(dgu_summary) <- NULL
  return(dgu_summary)
}

