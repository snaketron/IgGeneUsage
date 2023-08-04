
# Description:
# compute pmax -> probability of DGU
get_pmax_dgu <- function(glm_ext) {
  
  get_pmax_gene <- function(x, beta_data) {
    p <- beta_data[,x]
    l <- base::length(p)
    o <- base::max(base::sum(p < 0)/l, base::sum(p>0)/l)
    return(o)
  }
  
  beta_data <- glm_ext$beta_gene_mu
  pmax <- base::vapply(X = base::seq_len(length.out = base::ncol(beta_data)),
                       FUN = get_pmax_gene,
                       beta_data = beta_data,
                       FUN.VALUE = base::numeric(length = 1))
  pmax <- 2 * pmax - 1
  return(pmax)
}

# Description:
# compute pmax -> from numeric vector
get_pmax <- function(x) {
  p <- x
  l <- base::length(p)
  o <- base::max(base::sum(p < 0)/l, base::sum(p>0)/l)
  return(o)
}
