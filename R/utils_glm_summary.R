
get_glm_summary_gu_anova <- function(glm, hdi_lvl, ud) {
  
  glm_summary <- rstan::summary(object = glm, 
                                digits = 4,
                                pars = "beta_gene_mu",
                                prob = c(0.5, (1-hdi_lvl)/2,
                                         1-(1-hdi_lvl)/2))
  
  glm_summary <- glm_summary$summary
  glm_summary <- base::data.frame(glm_summary)
  base::colnames(glm_summary) <- c("es_mean", "es_mean_se",
                                   "es_sd", "es_median",
                                   "es_L", "es_H",
                                   "Neff", "Rhat")
  glm_summary[, c("Rhat", "Neff")] <- NULL
  
  par <- base::rownames(glm_summary)
  par <- base::gsub(pattern = "beta_gene_mu|\\[|\\]", replacement = '', x = par)
  par <- base::do.call(rbind, base::strsplit(x = par, split = ','))
  
  glm_summary$gene_id <- base::as.numeric(par[,2])
  glm_summary$group_id <- base::as.numeric(par[,1])
  glm_summary$contrast <- ud$contrast
  
  # group map
  group_map <- base::data.frame(group_name = ud$group_names,
                          group_id = ud$group_id)
  group_map <- group_map[base::duplicated(group_map)==F,]
  base::rownames(group_map) <- group_map$group_id
  
  glm_summary$gene_name <- ud$gene_names[glm_summary$gene_id]
  glm_summary$condition <- group_map[base::as.character(
    glm_summary$group_id), "group_name"]
  glm_summary$pmax <- NA
  
  
  # get pmax
  glm_ext <- rstan::extract(object = glm, par = "beta_gene_mu")$beta_gene_mu
  for(i in base::seq_len(base::nrow(glm_summary))) {
    glm_summary$pmax[i] <- get_pmax(
      x = glm_ext[,glm_summary$group_id[i], glm_summary$gene_id[i]])
  }
  
  # remove unused vars
  glm_summary$gene_id <- NULL
  glm_summary$group_id <- NULL
  base::rownames(glm_summary) <- NULL
  return(glm_summary)
}


get_glm_summary_gu_univar <- function(glm, hdi_lvl, ud) {
  
  glm_summary <- rstan::summary(object = glm, 
                                digits = 4,
                                pars = "prob_gene",
                                prob = c(0.5, (1-hdi_lvl)/2,
                                         1-(1-hdi_lvl)/2))
  
  glm_summary <- glm_summary$summary
  glm_summary <- base::data.frame(glm_summary)
  base::colnames(glm_summary) <- c("prob_mean", "prob_mean_se",
                                   "prob_sd", "prob_median",
                                   "prob_L", "prob_H",
                                   "Neff", "Rhat")
  glm_summary[, c("Rhat", "Neff")] <- NULL
  
  par <- base::rownames(glm_summary)
  par <- base::gsub(pattern = "prob_gene|\\[|\\]", replacement = '', x = par)
  par <- base::do.call(rbind, base::strsplit(x = par, split = ','))
  
  glm_summary$gene_id <- base::as.numeric(par)
  glm_summary$gene_name <- ud$gene_names[glm_summary$gene_id]
  glm_summary$condition <- ud$group_names[1]
  
  # remove unused vars
  glm_summary$gene_id <- NULL
  base::rownames(glm_summary) <- NULL
  return(glm_summary)
}


get_glm_summary_dgu <- function(glm, hdi_lvl, ud) {
  glm_summary <- rstan::summary(object = glm, 
                                digits = 4,
                                pars = "beta_gene_mu",
                                prob = c(0.5, (1-hdi_lvl)/2,
                                         1-(1-hdi_lvl)/2))
  
  glm_summary <- glm_summary$summary
  glm_summary <- base::data.frame(glm_summary)
  base::colnames(glm_summary) <- c("es_mean", "es_mean_se",
                                   "es_sd", "es_median",
                                   "es_L", "es_H",
                                   "Neff", "Rhat")
  glm_summary[, c("Rhat", "Neff")] <- NULL
  glm_summary$contrast <- ud$contrast
  
  # extract data and compute pmax
  glm_ext <- rstan::extract(object = glm, par = "beta_gene_mu")
  glm_summary$pmax <- get_pmax_dgu(glm_ext = glm_ext)
  
  # add gene id
  glm_summary$gene_name <- ud$gene_names
  
  return(glm_summary)
}
