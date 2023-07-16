
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


# Wilcoxon U test
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

