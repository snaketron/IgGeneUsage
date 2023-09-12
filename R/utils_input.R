
# Description:
# Provided the input arguments, this function checks their
# validity. It stops the execution if a problem is encountered
# and prints out warnings.
check_dgu_input <- function(ud,
                            mcmc_chains,
                            mcmc_cores,
                            mcmc_steps,
                            mcmc_warmup,
                            hdi_lvl) {
  
  if(base::missing(ud) || base::is.null(ud) ||
     base::missing(mcmc_chains) || base::is.null(mcmc_chains) ||
     base::missing(mcmc_steps) || base::is.null(mcmc_steps) ||
     base::missing(mcmc_warmup) || base::is.null(mcmc_warmup) ||
     base::missing(mcmc_cores) || base::is.null(mcmc_cores) ||
     base::missing(hdi_lvl) || base::is.null(hdi_lvl)) {
    stop("arguments must be specified")
  }
  
  check_ud(ud = ud)
  check_mcmc_steps(mcmc_steps = mcmc_steps,
                   mcmc_warmup = mcmc_warmup)
  check_mcmc_chains(mcmc_chains = mcmc_chains)
  check_mcmc_cores(mcmc_cores = mcmc_cores)
  check_hdi(hdi_lvl = hdi_lvl)
}


# Description:
# MCMC Iterations check
check_mcmc_steps <- function(mcmc_steps,
                             mcmc_warmup) {
  
  if(base::length(mcmc_steps) != 1 | 
     base::length(mcmc_warmup) != 1) {
    stop("mcmc_steps >= 500 & mcmc_warmup >= 100.")
  }
  
  if(!base::is.numeric(mcmc_steps) | 
     !base::is.numeric(mcmc_warmup)) {
    stop("mcmc_steps >= 500 & mcmc_warmup >= 100.")
  }
  
  
  if(base::is.finite(x = mcmc_steps)==FALSE | 
     base::is.finite(x = mcmc_warmup)==FALSE) {
    stop("mcmc_steps >= 500 & mcmc_warmup >= 100.")
  }
  
  
  if(base::as.integer(x = mcmc_steps) < 500 | 
     base::as.integer(x = mcmc_warmup) < 100) {
    stop("mcmc_steps >= 500 & mcmc_warmup >= 100.")
  }
  
  
  if(base::as.integer(x = mcmc_steps) <= base::as.integer(x = mcmc_warmup)) {
    stop("mcmc_steps > mcmc_warmup")
  }
}

# Description:
# MCMC Chain number check
check_mcmc_chains <- function(mcmc_chains) {
  if(base::length(mcmc_chains) != 1) {
    stop("mcmc_chains must be a positive integer > 0")
  }
  
  if(!base::is.numeric(mcmc_chains)) {
    stop("mcmc_chains must be a positive integer > 0")
  }
  
  if(base::is.finite(x = mcmc_chains) == FALSE) {
    stop("mcmc_chains must be a positive integer > 0")
  }
  
  if(base::as.integer(x = mcmc_chains) <= 0) {
    stop("mcmc_chains must be a positive integer > 0")
  }
}

# Description:
# MCMC Cores number check
check_mcmc_cores <- function(mcmc_cores) {
  if(base::length(mcmc_cores) != 1) {
    stop("mcmc_cores must be a positive integer > 0")
  }
  
  if(base::is.numeric(mcmc_cores) == FALSE) {
    stop("mcmc_cores must be a positive integer > 0")
  }
  
  if(base::is.finite(x = mcmc_cores) == FALSE) {
    stop("mcmc_cores must be a positive integer > 0")
  }
  
  if(base::as.integer(x = mcmc_cores) <= 0) {
    stop("mcmc_cores must be a positive integer > 0")
  }
}

# Description:
# HDI input check
check_hdi <- function(hdi_lvl) {
  if(base::length(hdi_lvl) != 1) {
    stop('hdi_lvl must be a number in range (0, 1)')
  }
  
  if(base::is.numeric(hdi_lvl) == FALSE) {
    stop('hdi_lvl must be a number in range (0, 1)')
  }
  
  if(hdi_lvl >= 1 | hdi_lvl <= 0) {
    stop('hdi_lvl must be a number in range (0, 1)')
  }
}

# Description:
# checks ud data.frame
check_ud <- function(ud) {
  
  if(base::is.data.frame(ud) == FALSE) {
    stop("ud must be data.frame")
  }
  
  if(base::ncol(ud) != 4) {
    stop("ud must contain the following columns: 'sample_id',
         'condition', 'gene_name' and 'gene_usage_count'")
  }
  
  cols <- c("sample_id", "condition", "gene_name", "gene_usage_count")
  if(base::all(base::colnames(ud) %in% cols) == FALSE) {
    stop("ud must contain the following columns: 'sample_id',
         'condition', 'gene_name' and 'gene_usage_count'")
  }
  
  if(base::nrow(ud) <= 1) {
    stop("ud must contain at least 2 data points")
  }
  
  if(base::typeof(ud$sample_id) != "character") {
    stop("column sample_id must be of character type.")
  }
  
  if(base::typeof(ud$condition) != "character") {
    stop("column condition must be of character type.")
  }
  
  if(base::typeof(ud$gene_name) != "character") {
    stop("column gene_name must be of character type.")
  }
  
  if(base::typeof(ud$gene_usage_count) != "numeric" &
     base::typeof(ud$gene_usage_count) != "double" &
     base::typeof(ud$gene_usage_count) != "integer") {
    stop("column gene_usage_count must be of numeric type.")
  }
}
