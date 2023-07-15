# Description:
# From the input data format guess the type of the analysis: paired or unpaired
get_analysis_type <- function(u) {
  if(all(c("gene_usage_count_1", "gene_usage_count_2") %in% colnames(u))) {
    return("paired")
  }
  return("unpaired")
}


# Description:
# Get usage data
get_usage <- function(u) {
  analysis_type <- get_analysis_type(u)
  if(analysis_type=="unpaired") {
    return(get_unpaired_usage(u = u))
  }
  if(analysis_type=="paired") {
    return(get_paired_usage(u = u))
  }
}


# Description:
# Parse input data for unpaired analysis
get_unpaired_usage <- function(u) {
  # if the same sample_id is present in both conditions
  k <- base::paste0(u$sample_id, '_', u$condition)
  if(base::length(base::unique(k)) != base::length(base::unique(u$sample_id))){
    warning("Same sample_id in both conditions, sample_id's extra coded")
    u$sample_id <- k
    rm(k)
  }
  
  # makes sure that we add usage = 0 for missing genes of certain samples 
  u <- tidyr::complete(u, sample_id, gene_name, 
                       fill = base::list(gene_usage_count = 0))
  
  # format data
  n <- stats::aggregate(gene_usage_count~sample_id+condition, 
                        FUN = sum, data = u)
  n$total_usage_count <- n$gene_usage_count
  n$gene_usage_count <- NULL
  
  u <- merge(x = u, y = n, by = c("sample_id", "condition"))
  u$gene_usage_prop <- u$gene_usage_count/u$total_usage_count
  rm(n)
  
  # get Y matrix
  Y <- reshape2::acast(data = u, 
                       formula = gene_name~sample_id,
                       drop = FALSE, 
                       value.var = "gene_usage_count",
                       fill = 0, 
                       fun.aggregate = sum)
  sample_ids <- base::colnames(Y)
  gene_names <- base::rownames(Y)
  
  N <- base::apply(X = Y, MARGIN = 2, FUN = sum)
  N <- base::as.numeric(N)
  
  X <- base::numeric(length = length(N))
  
  # sorted unique conditions
  cs <- base::sort(unique(u$condition)[1], decreasing = TRUE)
  # design vector X
  X <- base::ifelse(test = sample_ids %in% u$sample_id[u$codition == cs[1]], 
                    yes = +1, no = -1)
  X_org <- ifelse(test = X==1, yes = cs[1], no = cs[2])
  
  # add contrast
  contrast <- paste0(unique(X_org[X == +1]), " - ", unique(X_org[X == -1]))
  
  return (base::list(Y = Y, 
                     N = base::as.numeric(N), 
                     N_sample = base::ncol(Y), 
                     N_gene = base::nrow(Y),
                     X = X, 
                     X_org = X_org, 
                     gene_names = gene_names,
                     sample_names = sample_ids, 
                     proc_ud = u,
                     contrast = ))
}



# Description:
# Parse input data for paired analysis
get_paired_usage <- function(u) {
  
  # makes sure that we add usage = 0 for missing genes of certain samples 
  u <- tidyr::complete(u, sample_id, gene_name, 
                       fill = base::list(gene_usage_count_1 = 0,
                                         gene_usage_count_2 = 0))
  
  # format data
  u_1 <- aggregate(gene_usage_count_1~sample_id, FUN = sum, data = u)
  u_1$total_usage_count_1 <- u_1$gene_usage_count_1
  u_1$gene_usage_count_1 <- NULL
  
  u_2 <- aggregate(gene_usage_count_2~sample_id, FUN = sum, data = u)
  u_2$total_usage_count_2 <- u_2$gene_usage_count_2
  u_2$gene_usage_count_2 <- NULL
  
  u <- merge(x = u, y = u_1, by = c("sample_id"))
  u <- merge(x = u, y = u_2, by = c("sample_id"))
  u$gene_usage_prop_1 <- u$gene_usage_count_1/u$total_usage_count_1
  u$gene_usage_prop_2 <- u$gene_usage_count_2/u$total_usage_count_2
  rm(u_1, u_2)
  
  # get matrices
  # Y_1
  Y_1 <- reshape2::acast(data = u, 
                         formula = gene_name~sample_id,
                         drop = FALSE, 
                         value.var = "gene_usage_count_1",
                         fill = 0, 
                         fun.aggregate = base::sum)
  Y_1 <- Y_1[rownames(Y_1),]
  
  
  # Y_2
  Y_2 <- reshape2::acast(data = u, 
                         formula = gene_name~sample_id,
                         drop = FALSE, 
                         value.var = "gene_usage_count_2",
                         fill = 0, 
                         fun.aggregate = base::sum)
  # order Y_1 and Y_2 by same gene names
  Y_2 <- Y_2[rownames(Y_1),] 
  
  
  # compute total usage
  N_1 <- apply(X = Y_1, MARGIN = 2, FUN = sum)
  N_2 <- apply(X = Y_2, MARGIN = 2, FUN = sum)
  
  
  return (base::list(Y_1 = Y_1, 
                     Y_2 = Y_2, 
                     N_1 = N_1,
                     N_2 = N_2,
                     N_sample = base::ncol(Y_1), 
                     N_gene = base::nrow(Y_1), 
                     gene_names = base::rownames(Y_1),
                     sample_names = base::colnames(Y_1), 
                     proc_ud = u,
                     contrast = NA))
}


# Description:
# compute pmax -> probability of DGU
get_pmax <- function(glm_ext) {
  
  getPmaxGene <- function(x, beta_data) {
    p <- beta_data[,x]
    l <- length(p)
    o <- max(sum(p < 0)/l, sum(p>0)/l)
    return(o)
  }
  
  beta_data <- glm_ext$beta_gene_mu
  pmax <- vapply(X = seq_len(length.out = ncol(beta_data)),
                 FUN = getPmaxGene,
                 beta_data = beta_data,
                 FUN.VALUE = numeric(length = 1))
  pmax <- 2 * pmax - 1
  return(pmax)
}


# Description
# Posterior-predictive check in repertoires
get_ppc_rep <- function(glm,
                        ud,
                        hdi_lvl) {
  
  # summaries
  yhat <- summary(object = glm, pars = c("Yhat", "Yhat_rep"),
                  prob = c(0.5, (1-hdi_lvl)/2, 1-(1-hdi_lvl)/2))
  yhat <- yhat$summary
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
  yhat_c <- yhat[yhat$par_name == "Yhat", ]
  yhat_p <- yhat[yhat$par_name == "Yhat_rep", ]
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
    yhat$observed_count[i] <- ud$Y[yhat$G[i], yhat$R[i]]
    yhat$observed_prop[i] <- yhat$observed_count[i]/ud$N[yhat$R[i]]
    yhat$condition[i] <- ud$Xorg[yhat$sample_name[i]]
  }
  return (yhat)
}


# Description
# Posterior-predictive in a biological condition
get_ppc_condition <- function(glm,
                              ud,
                              hdi_lvl) {
  
  # summaries
  yhat <- summary(object = glm, 
                  pars = c("Yhat_condition"),
                  prob = c(0.5, (1-hdi_lvl)/2, 1-(1-hdi_lvl)/2))
  yhat <- yhat$summary
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
  yhat$X <- par_i[, 1]
  yhat$G <- par_i[, 2]
  rm(par_i)
  
  yhat$par <- NULL
  yhat$par_name <- NULL
  colnames(yhat)[base::seq_len(length.out = 6)] <- paste(
    "ppc", colnames(yhat)[base::seq_len(length.out = 6)], "prop", sep = "_")
  
  yhat$X <- ifelse(test = yhat$X == 1, yes = 1, no = -1)
  yhat$condition <- NA
  yhat$gene_name <- NA
  yhat$observed_prop <- NA
  for(i in base::seq_len(length.out = nrow(yhat))) {
    yhat$gene_name[i] <- ud$gene_names[yhat$G[i]]
    yhat$observed_prop[i] <- mean(ud$Y[yhat$G[i], ]/ud$N)
    yhat$condition[i] <- ud$Xorg[ud$X == yhat$X[i]][1]
  }
  
  return (yhat)
}


# t.test
get_ttest <- function(ud, 
                      paired) {
  
  get_ttest_run <- function(x, Ys, Xs, Ns, paired) {
    return(try(stats::t.test((Ys[x, ]/Ns)~X, 
                             paired = paired), 
               silent = TRUE))
  }
  
  get_ttest_summary <- function(x) {
    if(inherits(x = x, what = 'try-error') == TRUE) {
      return(data.frame(t_test_pvalue = NA,
                        t_test_tvalue = NA,
                        t_test_L95 = NA,
                        t_test_H95 = NA,
                        stringsAsFactors = FALSE))
    }
    return(data.frame(t_test_pvalue = x$p.value,
                      t_test_tvalue = x$statistic,
                      t_test_L95 = x$conf.int[1],
                      t_test_H95 = x$conf.int[2],
                      stringsAsFactors = FALSE))
  }
  
  tout <- lapply(X = seq_len(length.out = ud$N_gene),
                 FUN = get_ttest_run,
                 Ys = ud$Y,
                 Xs = ud$X,
                 Ns = ud$N,
                 paired = paired)
  
  tout.summary <- do.call(rbind, lapply(tout, get_ttest_summary))
  tout.summary$gene_name <- ud$gene_names
  
  # multiple correction
  tout.summary$t_test_fdr_pvalue <- stats::p.adjust(
    p = tout.summary$t_test_pvalue, method = "fdr")
  
  return (tout.summary)
}

#
get_manu <- function(ud, 
                     paired) {
  
  get_manu_run <- function(x, Ys, Xs, Ns, paired) {
    return(try(stats::wilcox.test((
      Ys[x, ]/Ns)~Xs, paired = paired), silent = TRUE))
  }
  
  get_manu_summary <- function(x) {
    if(inherits(x = x, what = "try-error") == TRUE) {
      return(data.frame(u_test_pvalue = NA,
                        u_test_wvalue = NA,
                        stringsAsFactors = FALSE))
    }
    return(data.frame(u_test_pvalue = x$p.value,
                      u_test_wvalue = x$statistic,
                      stringsAsFactors = FALSE))
  }
  
  mout <- lapply(X = seq_len(length.out = ud$N_gene),
                 FUN = get_manu_run,
                 Ys = ud$Y,
                 Xs = ud$X,
                 Ns = ud$N,
                 paired = paired)
  
  mout_summary <- do.call(rbind, lapply(X = mout, FUN = get_manu_summary))
  mout_summary$gene_name <- ud$gene_names
  
  # multiple correction
  mout_summary$u_test_fdr_pvalue <- stats::p.adjust(
    p = mout_summary$u_test_pvalue, method = "fdr")
  
  return (mout_summary)
}

