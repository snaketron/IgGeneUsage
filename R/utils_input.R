
# Description:
# Provided the input arguments, this function checks their
# validity. It stops the execution if a problem is encountered
# and prints out warnings.
check_dgu_input <- function(ud,
                            mcmc_chains,
                            mcmc_cores,
                            mcmc_steps,
                            mcmc_warmup,
                            hdi_lvl,
                            paired) {
  
  if(missing(ud) || is.null(ud) ||
     missing(mcmc_chains) || is.null(mcmc_chains) ||
     missing(mcmc_steps) || is.null(mcmc_steps) ||
     missing(mcmc_warmup) || is.null(mcmc_warmup) ||
     missing(mcmc_cores) || is.null(mcmc_cores) ||
     missing(hdi_lvl) || is.null(hdi_lvl) ||
     missing(paired) || is.null(paired)) {
    stop("arguments must be specified")
  }
  
  check_ud(ud = ud)
  check_mcmc_steps(mcmc_steps = mcmc_steps,
                   mcmc_warmup = mcmc_warmup)
  check_mcmc_chains(mcmc_chains = mcmc_chains)
  check_mcmc_cores(mcmc_cores = mcmc_cores)
  check_hdi(hdi_lvl = hdi_lvl)
  check_paired(paired = paired)
}


# Description:
# MCMC Iterations check
check_mcmc_steps <- function(mcmc_steps,
                             mcmc_warmup) {
  
  if(length(mcmc_steps) != 1 | 
     length(mcmc_warmup) != 1) {
    stop("mcmc_steps >= 500 & mcmc_warmup >= 100.")
  }
  
  if(!is.numeric(mcmc_steps) | 
     !is.numeric(mcmc_warmup)) {
    stop("mcmc_steps >= 500 & mcmc_warmup >= 100.")
  }
  
  
  if(is.finite(x = mcmc_steps)==FALSE | 
     is.finite(x = mcmc_warmup)==FALSE) {
    stop("mcmc_steps >= 500 & mcmc_warmup >= 100.")
  }
  
  
  if(as.integer(x = mcmc_steps) < 500 | 
     as.integer(x = mcmc_warmup) < 100) {
    stop("mcmc_steps >= 500 & mcmc_warmup >= 100.")
  }
  
  
  if(as.integer(x = mcmc_steps) <= as.integer(x = mcmc_warmup)) {
    stop("mcmc_steps > mcmc_warmup")
  }
}

# Description:
# MCMC Chain number check
check_mcmc_chains <- function(mcmc_chains) {
  if(length(mcmc_chains) != 1) {
    stop("mcmc_chains must be a positive integer > 0")
  }
  
  if(!is.numeric(mcmc_chains)) {
    stop("mcmc_chains must be a positive integer > 0")
  }
  
  if(is.finite(x = mcmc_chains) == FALSE) {
    stop("mcmc_chains must be a positive integer > 0")
  }
  
  if(as.integer(x = mcmc_chains) <= 0) {
    stop("mcmc_chains must be a positive integer > 0")
  }
}

# Description:
# MCMC Cores number check
check_mcmc_cores <- function(mcmc_cores) {
  if(length(mcmc_cores) != 1) {
    stop("mcmc_cores must be a positive integer > 0")
  }
  
  if(is.numeric(mcmc_cores) == FALSE) {
    stop("mcmc_cores must be a positive integer > 0")
  }
  
  if(is.finite(x = mcmc_cores) == FALSE) {
    stop("mcmc_cores must be a positive integer > 0")
  }
  
  if(as.integer(x = mcmc_cores) <= 0) {
    stop("mcmc_cores must be a positive integer > 0")
  }
}

# Description:
# HDI input check
check_hdi <- function(hdi_lvl) {
  if(length(hdi_lvl) != 1) {
    stop('hdi_lvl must be a number in range (0, 1)')
  }
  
  if(is.numeric(hdi_lvl) == FALSE) {
    stop('hdi_lvl must be a number in range (0, 1)')
  }
  
  if(hdi_lvl >= 1 | hdi_lvl <= 0) {
    stop('hdi_lvl must be a number in range (0, 1)')
  }
}

# Description:
# checks ud data.frame
check_ud <- function(ud) {
  
  if(is.data.frame(ud) == FALSE) {
    stop("ud must be data.frame")
  }
  
  cols <- c("individual_id", "condition", "gene_name", "gene_usage_count")
  cols_r <- c("individual_id", "condition", "gene_name", "gene_usage_count", 
              "replicate")
  if(all(colnames(ud) %in% cols | colnames(ud) %in% cols_r) == FALSE) {
    stop("ud must contain the following columns: 'individual_id',
         'condition', 'gene_name' and 'gene_usage_count' and optionally 
         'replicate'")
  }
  
  if(nrow(ud) <= 1) {
    stop("ud must contain at least 2 data points")
  }
  
  if(typeof(ud$individual_id) != "character") {
    stop("column individual_id must be character")
  }
  
  if(typeof(ud$condition) != "character") {
    stop("column condition must be character")
  }
  
  if(typeof(ud$gene_name) != "character") {
    stop("column gene_name must be character")
  }
  
  if(typeof(ud$gene_usage_count) != "numeric" &
     typeof(ud$gene_usage_count) != "double" &
     typeof(ud$gene_usage_count) != "integer") {
    stop("column gene_usage_count must be numeric")
  }
  
  if("replicate" %in% colnames(ud)) {
    if(typeof(ud$replicate) != "character" & typeof(ud$replicate) != "numeric"){
      stop("column replicate must be character or numeric")
    }
    
    k <- ud[duplicated(ud[,c("individual_id","condition","replicate")])==FALSE,
            c("individual_id","condition")]
    if(all(table(apply(X = k, MARGIN = 1, FUN = paste0, collapse = '|'))==1)){
      warning("one replicate available, no-replicates model used")
    }
  }
}



# Description:
# paired input check
check_paired <- function(paired) {
  if(length(paired) != 1) {
    stop('paired must be a logical (TRUE or FALSE)')
  }
  
  if(is.logical(paired) == FALSE) {
    stop('paired must be a logical (TRUE or FALSE)')
  }
}