

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


  if(missing(usage.data) || is.null(usage.data) ||
     missing(mcmc.chains) || is.null(mcmc.chains) ||
     missing(mcmc.steps) || is.null(mcmc.steps) ||
     missing(mcmc.warmup) || is.null(mcmc.warmup) ||
     missing(mcmc.cores) || is.null(mcmc.cores) ||
     missing(hdi.level) || is.null(hdi.level)) {
    stop("arguments must be specified")
  }


  checkUsageData(usage.data = usage.data)
  checkMcmcSteps(mcmc.steps = mcmc.steps,
                 mcmc.warmup = mcmc.warmup)
  checkMcmcChains(mcmc.chains = mcmc.chains)
  checkMcmcCores(mcmc.cores = mcmc.cores)
  checkHdi(hdi.level = hdi.level)
}


# Description:
# Usage data check if data.frame
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
  
  if(length(unique(usage.data$condition)) != 2) {
    stop("exactly 2 biological conditions must be provided.")
  }
}






# Description:
# MCMC Iterations check
checkMcmcSteps <- function(mcmc.steps,
                           mcmc.warmup) {

  if(length(mcmc.steps) != 1 | length(mcmc.warmup) != 1) {
    stop("mcmc.steps >= 500 & mcmc.warmup >= 100.")
  }

  if(!is.numeric(mcmc.steps) | !is.numeric(mcmc.warmup)) {
    stop("mcmc.steps >= 500 & mcmc.warmup >= 100.")
  }


  if(is.finite(x = mcmc.steps)==FALSE | is.finite(x = mcmc.warmup)==FALSE) {
    stop("mcmc.steps >= 500 & mcmc.warmup >= 100.")
  }


  if(as.integer(x = mcmc.steps) < 500 | as.integer(x = mcmc.warmup) < 100) {
    stop("mcmc.steps >= 500 & mcmc.warmup >= 100.")
  }


  if(as.integer(x = mcmc.steps) <= as.integer(x = mcmc.warmup)) {
    stop("mcmc.steps > mcmc.warmup")
  }
}


# Description:
# MCMC Chain number check
checkMcmcChains <- function(mcmc.chains) {
  # CHECK: mcmc.chains
  if(length(mcmc.chains) != 1) {
    stop("mcmc.chains must be a positive integer > 0")
  }

  if(!is.numeric(mcmc.chains)) {
    stop("mcmc.chains must be a positive integer > 0")
  }

  if(is.finite(x = mcmc.chains) == FALSE) {
    stop("mcmc.chains must be a positive integer > 0")
  }

  if(as.integer(x = mcmc.chains) <= 0) {
    stop("mcmc.chains must be a positive integer > 0")
  }
}


# Description:
# MCMC Cores number check
checkMcmcCores <- function(mcmc.cores) {
  if(length(mcmc.cores) != 1) {
    stop("mcmc.cores must be a positive integer > 0")
  }

  if(is.numeric(mcmc.cores) == FALSE) {
    stop("mcmc.cores must be a positive integer > 0")
  }

  if(is.finite(x = mcmc.cores) == FALSE) {
    stop("mcmc.cores must be a positive integer > 0")
  }

  if(as.integer(x = mcmc.cores) <= 0) {
    stop("mcmc.cores must be a positive integer > 0")
  }
}


# Description:
# HDI input check
checkHdi <- function(hdi.level) {
  if(length(hdi.level) != 1) {
    stop('hdi.level must be a number in range (0, 1)')
  }

  if(is.numeric(hdi.level) == FALSE) {
    stop('hdi.level must be a number in range (0, 1)')
  }

  if(hdi.level >= 1 | hdi.level <= 0) {
    stop('hdi.level must be a number in range (0, 1)')
  }
}


