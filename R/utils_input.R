

# Description:
# Provided the input arguments, this function checks their
# validity. It stops the execution if a problem is encountered
# and prints out warnings.
check_input <- function(ud,
                        mcmc_chains,
                        mcmc_cores,
                        mcmc_steps,
                        mcmc_warmup,
                        hdi_lvl) {
  
  
  if(missing(ud) || is.null(ud) ||
     missing(mcmc_chains) || is.null(mcmc_chains) ||
     missing(mcmc_steps) || is.null(mcmc_steps) ||
     missing(mcmc_warmup) || is.null(mcmc_warmup) ||
     missing(mcmc_cores) || is.null(mcmc_cores) ||
     missing(hdi_lvl) || is.null(hdi_lvl)) {
    stop("arguments must be specified")
  }
  
  analysis_type <- get_analysis_type(ud)
  if(analysis_type == "paired") {
    check_usage_data_paired(ud = ud)
  }
  if(analysis_type == "unpaired") {
    check_usage_data_unpaired(ud = ud)
  }
  check_mcmc_steps(mcmc_steps = mcmc_steps,
                   mcmc_warmup = mcmc_warmup)
  check_mcmc_chains(mcmc_chains = mcmc_chains)
  check_mcmc_cores(mcmc_cores = mcmc_cores)
  check_hdi(hdi_lvl = hdi_lvl)
}


# Description:
# Provided the input arguments, this function checks their
# validity. It stops the execution if a problem is encountered
# and prints out warnings.
check_gu_input <- function(ud,
                           mcmc_chains,
                           mcmc_cores,
                           mcmc_steps,
                           mcmc_warmup,
                           hdi_lvl) {
  
  
  if(missing(ud) || is.null(ud) ||
     missing(mcmc_chains) || is.null(mcmc_chains) ||
     missing(mcmc_steps) || is.null(mcmc_steps) ||
     missing(mcmc_warmup) || is.null(mcmc_warmup) ||
     missing(mcmc_cores) || is.null(mcmc_cores) ||
     missing(hdi_lvl) || is.null(hdi_lvl)) {
    stop("arguments must be specified")
  }
  
  analysis_type <- get_analysis_type(ud)
  if(analysis_type == "paired") {
    stop("Gene usage (GU) analysis requires table with 4 columns: 
    'sample_id', 'gene_name', 'condition', 'gene_usage_count'")
  }
  if(analysis_type == "unpaired") {
    check_usage_data_unpaired(ud = ud)
  }
  check_mcmc_steps(mcmc_steps = mcmc_steps,
                   mcmc_warmup = mcmc_warmup)
  check_mcmc_chains(mcmc_chains = mcmc_chains)
  check_mcmc_cores(mcmc_cores = mcmc_cores)
  check_hdi(hdi_lvl = hdi_lvl)
}


# Description:
# Usage data check if data.frame
check_usage_data_unpaired <- function(ud) {
  
  if(is.data.frame(ud) == FALSE) {
    stop("ud must be data.frame")
  }
  
  if(ncol(ud) != 4) {
    stop("ud must contain the following columns: 'sample_id',
         'condition', 'gene_name' and 'gene_usage_count'")
  }
  
  cols <- c("sample_id", "condition", "gene_name", "gene_usage_count")
  if(all(colnames(ud) %in% cols) == FALSE) {
    stop("ud must contain the following columns: 'sample_id',
         'condition', 'gene_name' and 'gene_usage_count'")
  }
  
  if(nrow(ud) <= 1) {
    stop("ud must contain at least 2 data points")
  }
  
  if(typeof(ud$sample_id) != "character") {
    stop("column sample_id must be of character type.")
  }
  
  if(typeof(ud$condition) != "character") {
    stop("column condition must be of character type.")
  }
  
  if(typeof(ud$gene_name) != "character") {
    stop("column gene_name must be of character type.")
  }
  
  if(typeof(ud$gene_usage_count) != "numeric" &
     typeof(ud$gene_usage_count) != "double" &
     typeof(ud$gene_usage_count) != "integer") {
    stop("column gene_usage_count must be of numeric type.")
  }
  
  if(length(unique(ud$condition)) != 2) {
    stop("exactly 2 biological conditions must be provided.")
  }
  
  if(min(length(unique(ud$sample_id[ud$condition==unique(
    ud$condition)[1]])),
    length(unique(ud$sample_id[ud$condition==unique(
      ud$condition)[2]])))==1) {
    warning("no replicates provided for at least one of the conditions.")
  }
}

# Description:
# Usage data check if data.frame
check_usage_data_paired <- function(ud) {
  
  if(is.data.frame(ud) == FALSE) {
    stop("ud must be data.frame")
  }
  
  if(ncol(ud) != 4) {
    stop("ud must contain the following columns: 'sample_id',
         'gene_name' and 'gene_usage_count_1' and gene_usage_count_2")
  }
  
  cols <- c("sample_id", "gene_name", 
            "gene_usage_count_1", 
            "gene_usage_count_2")
  if(all(colnames(ud) %in% cols) == FALSE) {
    stop("ud must contain the following columns: 'sample_id',
         'gene_name' and 'gene_usage_count_1' and gene_usage_count_2")
  }
  
  if(nrow(ud) <= 1) {
    stop("ud must contain at least 2 data points")
  }
  
  if(typeof(ud$sample_id) != "character") {
    stop("column sample_id must be of character type.")
  }
  
  if(typeof(ud$gene_name) != "character") {
    stop("column gene_name must be of character type.")
  }
  
  if(typeof(ud$gene_usage_count_1) != "numeric" &
     typeof(ud$gene_usage_count_1) != "double" &
     typeof(ud$gene_usage_count_1) != "integer") {
    stop("column gene_usage_count_1 must be of numeric type.")
  }
  
  if(typeof(ud$gene_usage_count_2) != "numeric" &
     typeof(ud$gene_usage_count_2) != "double" &
     typeof(ud$gene_usage_count_2) != "integer") {
    stop("column gene_usage_count_2 must be of numeric type.")
  }
}

# Description:
# MCMC Iterations check
check_mcmc_steps <- function(mcmc_steps,
                             mcmc_warmup) {
  
  if(length(mcmc_steps) != 1 | length(mcmc_warmup) != 1) {
    stop("mcmc_steps >= 500 & mcmc_warmup >= 100.")
  }
  
  if(!is.numeric(mcmc_steps) | !is.numeric(mcmc_warmup)) {
    stop("mcmc_steps >= 500 & mcmc_warmup >= 100.")
  }
  
  
  if(is.finite(x = mcmc_steps)==FALSE | is.finite(x = mcmc_warmup)==FALSE) {
    stop("mcmc_steps >= 500 & mcmc_warmup >= 100.")
  }
  
  
  if(as.integer(x = mcmc_steps) < 500 | as.integer(x = mcmc_warmup) < 100) {
    stop("mcmc_steps >= 500 & mcmc_warmup >= 100.")
  }
  
  
  if(as.integer(x = mcmc_steps) <= as.integer(x = mcmc_warmup)) {
    stop("mcmc_steps > mcmc_warmup")
  }
}


# Description:
# MCMC Chain number check
check_mcmc_chains <- function(mcmc_chains) {
  # CHECK: mcmc_chains
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


