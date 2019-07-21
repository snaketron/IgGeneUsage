

# Description:
# Provided the input arguments, this function checks their
# validity. It stops the execution if a problem is encountered
# and prints out warnings.
checkInput <- function(usage.data,
                       mcmc.chains,
                       mcmc.cores,
                       mcmc.steps,
                       mcmc.warmup,
                       hdi.level) {


  if(missing(usage.data) || is.null(usage.data) || is.na(usage.data) ||
     missing(mcmc.chains) || is.null(mcmc.chains) || is.na(mcmc.chains) |
     missing(mcmc.steps) || is.null(mcmc.steps) || is.na(mcmc.steps) ||
     missing(mcmc.warmup) || is.null(mcmc.warmup) || is.na(mcmc.warmup) |
     missing(mcmc.cores) || is.null(mcmc.cores) || is.na(mcmc.cores) ||
     missing(hdi.level) || is.null(hdi.level) || is.na(hdi.level)) {
    stop("arguments must be specified")
  }


  checkUsageData(usage.data = usage.data)
  checkMcmcIterations(mcmc.steps = mcmc.steps)
  checkMcmcWarmup(mcmc.warmup = mcmc.warmup)
  checkMcmcChains(mcmc.chains = mcmc.chains)
  checkMcmcCores(mcmc.cores = mcmc.cores)
  checkHdi(hdi.level = hdi.level)
}


# Description:
# Usage data check
checkUsageData <- function(usage.data) {

  if(is.data.frame(usage.data) == FALSE) {
    stop("usage.data must be data.frame")
  }

  if(ncol(usage.data) != 4) {
    stop("usage.data must contain the following columns: 'sample_id',
         'condition', 'gene_name' and 'gene_usage_count'")
  }

  correct.col <- c("sample_id", "condition",
                   "gene_name", "gene_usage_count")
  if(all(colnames(usage.data) %in% correct.col) == FALSE) {
    stop("usage.data must contain the following columns: 'sample_id',
         'condition', 'gene_name' and 'gene_usage_count'")
  }

  if(nrow(usage.data) <= 1) {
    stop("usage.data must contain at least 2 data points")
  }

  if(typeof(usage.data$sample_id) != "character") {
    stop("column sample_id must be of character type.")
  }

  if(typeof(usage.data$condition) != "character") {
    stop("column condition must be of character type.")
  }

  if(typeof(usage.data$gene_name) != "character") {
    stop("column gene_name must be of character type.")
  }

  if(typeof(usage.data$gene_usage_count) != "numeric" &
     typeof(usage.data$gene_usage_count) != "double" &
     typeof(usage.data$gene_usage_count) != "integer") {
    stop("column gene_usage_count must be of numeric type.")
  }
}


# Description:
# MCMC Iterations check
checkMcmcIterations <- function(mcmc.steps) {
  # CHECK: mcmc.steps
  if(length(mcmc.steps) != 1) {
    stop("the mcmc.steps must be a number > 0 (default = 10000).")
  }

  if(!is.numeric(mcmc.steps)) {
    stop("mcmc.steps must be a numeric argument (default = 10000).")
  }

  if(mcmc.steps <= 0) {
    stop("mcmc.steps must be larger than 0 (default = 10000).")
  }
}


# Description:
# MCMC Warmup check
checkMcmcWarmup <- function(mcmc.warmup) {
  # CHECK: mcmc.warmup
  if(length(mcmc.warmup) != 1) {
    stop("the mcmc.warmup must be a number > 0 (default = 5000).")
  }

  if(!is.numeric(mcmc.warmup)) {
    stop("mcmc.warmup must be a numeric argument (default = 5000).")
  }

  if(mcmc.warmup <= 0) {
    stop("mcmc.warmup must be larger than 0 (default = 5000).")
  }
}


# Description:
# MCMC Chain number check
checkMcmcChains <- function(mcmc.chains) {
  # CHECK: mcmc.chains
  if(length(mcmc.chains) != 1) {
    stop("mcmc.chains must be a positive integer > 0 (default = 1).")
  }

  if(!is.numeric(mcmc.chains)) {
    stop("mcmc.chains must be a positive integer > 0 (default = 1).")
  }

  if(mcmc.chains <= 0) {
    stop("mcmc.chains must be a positive integer > 0 (default = 1).")
  }
}


# Description:
# MCMC Cores number check
checkMcmcCores <- function(mcmc.cores) {
  # CHECK: cores
  if(length(mcmc.cores) != 1) {
    stop("mcmc.cores is numeric parameter.")
  }

  if(is.numeric(mcmc.cores) == FALSE) {
    stop("mcmc.cores is numeric parameter.")
  }

  if(mcmc.cores <= 0) {
    stop("mcmc.cores is numeric parameter >=1.")
  }
}


# Description:
# HDI input check
checkHdi <- function(hdi.level) {
  if(length(hdi.level) != 1) {
    stop("The HDI level must be in range (0, 1).")
  }

  if(is.numeric(hdi.level) == FALSE) {
    stop("The HDI level must be in range (0, 1).")
  }

  if(hdi.level >= 1 | hdi.level <= 0) {
    stop("The HDI level must be in range (0, 1).")
  }
}


