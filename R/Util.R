
# Description:
# Parse input data
get_usage <- function(u) {
  # if the same sample_id is present in both conditions
  key <- paste(u$sample_id, u$condition, sep = '_')
  if(length(unique(key)) != length(unique(u$sample_id))) {
    warning("Same sample_id in both conditions, sample_id's extra coded")
    u$sample_id <- key
    rm(key)
  }
  
  # get Y data, fill empty combinations with 0
  Y <- reshape2::acast(data = u, 
                       formula = gene_name~sample_id,
                       drop = FALSE, 
                       value.var = "gene_usage_count",
                       fill = 0, 
                       fun.aggregate = sum)
  
  sample_ids <- colnames(Y)
  gene_names <- rownames(Y)
  
  # get N data
  N_d <- stats::aggregate(gene_usage_count~sample_id, 
                          data = u,
                          FUN = sum, 
                          drop = FALSE)
  N <- N_d$gene_usage_count
  names(N) <- N_d$sample_id
  rm(N_d)
  N <- N[sample_ids]
  
  # get X data
  u <-  u[u$sample_id %in% sample_ids, ]
  u <- u[duplicated(u[, c("sample_id")]) == FALSE, ]
  
  X_d <- stats::aggregate(condition~sample_id, 
                          data = u, 
                          FUN = unique)
  X <- X_d$condition
  names(X) <- X_d$sample_id
  rm(X_d)
  X <- X[sample_ids]
  
  # get X mapping
  X_u <- sort(x = unique(X), decreasing = TRUE)
  x1 <- which(X == X_u[1])
  x2 <- which(X == X_u[2])
  X_m <- numeric(length = length(X))
  X_m[x1] <- 1
  X_m[x2] <- -1
  
  # compute processed usage data
  pu <- get_processed_ud(Y = Y, 
                         gene_names = gene_names, 
                         sample_ids = sample_ids, 
                         X = X)
  
  return (list(Y = Y, 
               N = N, 
               N_sample = ncol(Y), 
               N_gene = nrow(Y),
               X = X_m, 
               Xorg = X, 
               gene_names = gene_names,
               sample_names = sample_ids, 
               pu = pu))
}


# Description:
# Parse input data
get_paired_usage <- function(u) {
  u_t <- u[base::duplicated(u[, c("sample_id", "condition")])==FALSE,]
  u_t <- base::table(u_t$condition, u_t$sample_id)
  if(base::any(u_t!=1)) {
    stop(paste0("sample/s: ", 
                paste0(names(which(colSums(u_t)!=2)), collapse = ','),
                " not paired."))
  }
  
  # get Y data, fill empty combinations with 0
  Y <- reshape2::acast(data = u, 
                       formula = gene_name~sample_id+condition,
                       drop = FALSE, 
                       value.var = "gene_usage_count",
                       fill = 0, 
                       fun.aggregate = base::sum)
  Y <- reshape2::melt(Y)
  
  # processed data
  colnames(Y) <- c("gene_name", "sc", "gene_usage_count")
  k <- do.call(rbind, strsplit(x = as.character(Y$sc), split = '\\_'))
  Y$sample_id <- k[,1]
  Y$condition <- k[,2]
  Y$sc <- NULL
  
  # compute total usage
  N <- stats::aggregate(gene_usage_count~sample_id+condition, 
                        data = Y,
                        FUN = base::sum, 
                        drop = FALSE)
  N$total_usage_count <- N$gene_usage_count
  N$gene_usage_count <- NULL
  
  # merge usage and total usage
  Y <- base::merge(x = Y, y = N, 
                   by = c("sample_id", "condition"),
                   all.x = T)
  Y$gene_usage_prop <- Y$gene_usage_count/Y$total_usage_count
  Y$gene_name <- base::as.character(Y$gene_name)
  
  cs <- base::sort(base::unique(Y$condition))
  
  # get usage matrices
  Y_1 <- reshape2::acast(data = Y[Y$condition == cs[1], ], 
                         formula = gene_name~sample_id,
                         drop = FALSE, 
                         value.var = "gene_usage_count",
                         fill = 0, 
                         fun.aggregate = base::sum)
  Y_2 <- reshape2::acast(data = Y[Y$condition == cs[2], ], 
                         formula = gene_name~sample_id,
                         drop = FALSE, 
                         value.var = "gene_usage_count",
                         fill = 0, 
                         fun.aggregate = base::sum)
  Y_1 <- Y_1[base::sort(base::rownames(Y_1)), 
             base::sort(base::colnames(Y_1))]
  Y_2 <- Y_2[base::sort(base::rownames(Y_2)), 
             base::sort(base::colnames(Y_2))]
  
  N <- base::cbind(base::colSums(Y_1), 
                   base::colSums(Y_2))
  
  return (base::list(Y_1 = Y_1, 
                     Y_2 = Y_2, 
                     N = N, 
                     N_sample = base::nrow(Y), 
                     N_gene = base::nrow(Y_1), 
                     gene_names = base::rownames(Y_1),
                     sample_names = base::colnames(Y_1), 
                     pu = Y))
}


# Description:
# Takes an input ud of class SummarizedExperiment,
# and creates data.frame. Given that this function is used
# before the input control, it also does some rudimentary checks
get_SummarizedExperiment <- function(ud_se) {
  # get count data
  c_d <- SummarizedExperiment::assay(x = ud_se)
  c_d <- reshape2::melt(data = c_d, 
                        value.name = "gene_usage_count", 
                        as.is = TRUE)
  colnames(c_d) <- c("gene_name", "sample_id", "gene_usage_count")
  
  
  # get the sample (repertoire) data and convert it to data.frame
  coldata <- SummarizedExperiment::colData(x = ud_se)
  coldata <- base::as.data.frame(coldata)
  if(nrow(coldata) == 0 | ncol(coldata) == 0) {
    stop("colData(ud) is empty")
  }
  
  if(all(colnames(coldata) %in% c("condition", "sample_id")) == FALSE) {
    stop("colData(ud) must have 2 columns: 'condition' and 'sample_id'")
  }
  
  # merge both datasets into the standard 4-column IgGeneUsage input
  ud <- base::merge(x = c_d,
                    y = coldata[, c("condition", "sample_id")],
                    by = "sample_id", all = TRUE)
  
  return (ud)
}


get_pmax <- function(glm_ext) {
  
  getPmaxGene <- function(x, beta_data) {
    p <- beta_data[,x]
    l <- length(p)
    o <- max(sum(p < 0)/l, sum(p>0)/l)
    return(o)
  }
  
  beta_data <- glm_ext$beta_gene
  pmax <- vapply(X = seq_len(length.out = ncol(beta_data)),
                 FUN = getPmaxGene,
                 beta_data = beta_data,
                 FUN.VALUE = numeric(length = 1))
  
  # convert from [0.5, 1] -> [0, 1]
  pmax <- 2 * pmax - 1
  return(pmax)
}


# Description
# Posterior-predictive check within repertoires
get_ppc_rep <- function(glm,
                        ud,
                        hdi_lvl) {
  
  # summaries
  yhat <- summary(object = glm, pars = c("Yhat", "Yhat_individual"),
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
  yhat_p <- yhat[yhat$par_name == "Yhat_individual", ]
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
# Posterior-predictive check across genes
get_ppc_gene <- function(glm,
                         ud,
                         hdi_lvl) {
  
  # summaries
  yhat <- summary(object = glm, 
                  pars = c("Yhat_gene"),
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


# two sided t.test
get_ttest <- function(ud) {
  
  get_ttest_run <- function(x, Ys, Xs, Ns) {
    return(try(stats::t.test((Ys[x, ]/Ns)~Xs), silent = TRUE))
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
                 Ns = ud$N)
  
  tout.summary <- do.call(rbind, lapply(tout, get_ttest_summary))
  tout.summary$gene_name <- ud$gene_names
  
  # multiple correction
  tout.summary$t_test_fdr_pvalue <- stats::p.adjust(
    p = tout.summary$t_test_pvalue, method = "fdr")
  
  return (tout.summary)
}


get_manu <- function(ud) {
  
  get_manu_run <- function(x, Ys, Xs, Ns) {
    return(try(stats::wilcox.test((Ys[x, ]/Ns)~Xs), silent = TRUE))
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
                 Ns = ud$N)
  
  mout_summary <- do.call(rbind, lapply(X = mout, FUN = get_manu_summary))
  mout_summary$gene_name <- ud$gene_names
  
  # multiple correction
  mout_summary$u_test_fdr_pvalue <- stats::p.adjust(
    p = mout_summary$u_test_pvalue, method = "fdr")
  
  return (mout_summary)
}


get_processed_ud <- function(Y, 
                             gene_names, 
                             sample_ids, 
                             X) {
  # process usage
  pud <- vector(mode = "list", length = ncol(Y))
  for(i in base::seq_len(length.out = ncol(Y))) {
    pud[[i]] <- data.frame(gene_usage_count = Y[,i],
                           gene_name = gene_names,
                           total_usage_count = sum(Y[,i]),
                           gene_usage_prop = Y[,i]/sum(Y[,i]),
                           sample_id = sample_ids[i],
                           condition = X[i],
                           stringsAsFactors = FALSE)
  }
  pud <- do.call(rbind, pud)
  rownames(pud) <- NULL
  return(pud)
}

