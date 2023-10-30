get_contrast_map <- function(ud) {
  # group map
  group_map <- base::data.frame(group_name = ud$group_names,
                                group_id = ud$group_id)
  group_map <- group_map[base::duplicated(group_map)==FALSE,]
  base::rownames(group_map) <- group_map$group_id
  
  # contrast map
  c_n <- max(group_map$group_id)*(max(group_map$group_id)-1)/2
  c_map <- data.frame(c_id = base::seq_len(c_n))
  c_names <- c()
  for(i in 1:(max(group_map$group_id)-1)) {
    for(j in (i+1):max(group_map$group_id)) {
      c_names <- c(c_names,
                   paste0(group_map[base::as.character(i), "group_name"],
                          "-vs-",
                          group_map[base::as.character(j), "group_name"]))
    }
  }
  c_map$contrast <- c_names
  return(c_map)
}


get_dgu_summary <- function(glm, hdi_lvl, ud) {
  contrast_map <- get_contrast_map(ud)
  
  dgu_summary <- rstan::summary(object = glm, 
                                digits = 4,
                                pars = "dgu",
                                prob = c(0.5, (1-hdi_lvl)/2,
                                         1-(1-hdi_lvl)/2))
  
  dgu_summary <- dgu_summary$summary
  dgu_summary <- base::data.frame(dgu_summary)
  base::colnames(dgu_summary) <- c("es_mean", "es_mean_se",
                                   "es_sd", "es_median",
                                   "es_L", "es_H",
                                   "Neff", "Rhat")
  dgu_summary[, c("Rhat", "Neff")] <- NULL
  
  par <- base::rownames(dgu_summary)
  par <- base::gsub(pattern = "dgu|\\[|\\]", replacement = '', x = par)
  par <- base::do.call(rbind, base::strsplit(x = par, split = ','))
  
  
  dgu_summary$gene_id <- base::as.numeric(par[,1])
  dgu_summary$c_id <- base::as.numeric(par[,2])
  dgu_summary <- base::merge(x = dgu_summary, 
                             y = contrast_map, 
                             by = "c_id", 
                             all.x = TRUE)
  
  dgu_summary$gene_name <- ud$gene_names[dgu_summary$gene_id]
  dgu_summary$pmax <- NA
  
  # get pmax
  glm_ext <- rstan::extract(object = glm, par = "dgu")$dgu
  for(i in base::seq_len(base::nrow(dgu_summary))) {
    dgu_summary$pmax[i] <- get_pmax(
      x = glm_ext[,dgu_summary$gene_id[i], dgu_summary$c_id[i]])
  }
  
  # remove unused vars
  dgu_summary$gene_id <- NULL
  dgu_summary$c_id <- NULL
  base::rownames(dgu_summary) <- NULL
  return(dgu_summary)
}


get_dgu_prob_summary <- function(glm, hdi_lvl, ud) {
  contrast_map <- get_contrast_map(ud)
  
  dgu_summary <- rstan::summary(object = glm, 
                                digits = 4,
                                pars = "dgu_prob",
                                prob = c(0.5, (1-hdi_lvl)/2,
                                         1-(1-hdi_lvl)/2))
  
  dgu_summary <- dgu_summary$summary
  dgu_summary <- base::data.frame(dgu_summary)
  base::colnames(dgu_summary) <- c("es_prob_mean", "es_prob_mean_se",
                                   "es_prob_sd", "es_prob_median",
                                   "es_prob_L", "es_prob_H",
                                   "Neff", "Rhat")
  dgu_summary[, c("Rhat", "Neff")] <- NULL
  
  par <- base::rownames(dgu_summary)
  par <- base::gsub(pattern = "dgu_prob|\\[|\\]", replacement = '', x = par)
  par <- base::do.call(rbind, base::strsplit(x = par, split = ','))
  
  
  dgu_summary$gene_id <- base::as.numeric(par[,1])
  dgu_summary$c_id <- base::as.numeric(par[,2])
  dgu_summary <- base::merge(x = dgu_summary, 
                             y = contrast_map, 
                             by = "c_id", 
                             all.x = TRUE)
  
  dgu_summary$gene_name <- ud$gene_names[dgu_summary$gene_id]
  dgu_summary$pmax <- NA
  
  # get pmax
  glm_ext <- rstan::extract(object = glm, par = "dgu_prob")$dgu_prob
  for(i in base::seq_len(base::nrow(dgu_summary))) {
    dgu_summary$pmax[i] <- get_pmax(
      x = glm_ext[,dgu_summary$gene_id[i], dgu_summary$c_id[i]])
  }
  
  # remove unused vars
  dgu_summary$gene_id <- NULL
  dgu_summary$c_id <- NULL
  base::rownames(dgu_summary) <- NULL
  return(dgu_summary)
}

