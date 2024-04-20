
# Description
# Parse input data for GU analysis
get_usage <- function(u) {
  
  get_condition_of_sample <- function(sample_ids, 
                                      m) {
    # condition at sample
    condition_names <- character(length = length(sample_ids))
    for(i in 1:length(sample_ids)) {
      condition_names[i] <- m$condition[m$sample_id == sample_ids[i]]
    }
    condition_ids <- as.numeric(as.factor(condition_names))
    
    return(list(condition_ids = condition_ids,
                condition_names = condition_names))
  }
  
  get_condition_of_individual <- function(individual_ids, 
                                          individual_names, 
                                          m) {
    # condition at individual
    condition_names <- character(length = max(individual_ids))
    for(i in 1:max(individual_ids)) {
      q <- individual_names[individual_ids==i][1]
      condition_names[i] <- m$condition[m$individual_id == q][1]
    }
    condition_ids <- as.numeric(as.factor(condition_names))
    
    return(list(condition_ids = condition_ids,
                condition_names = condition_names))
  }
  
  check_paired <- function(u, 
                           has_balanced_replicates, 
                           has_replicates, 
                           has_conditions) {
    if(has_conditions==FALSE) {
      return(FALSE)
    }
    if(has_balanced_replicates==FALSE) {
      return(FALSE)
    }
    
    if(has_replicates) {
      return(has_balanced_replicates)
      # q <- u[duplicated(u[,c("individual_id","condition",
      #                        "replicate_id")]) == FALSE,]
      # q$f <- 1
      # q <- aggregate(f~individual_id+condition+replicate_id, 
      #                data = q, FUN = sum, drop = FALSE)
      # q$f[is.null(q$f)|is.na(q$f)] <- 0
      # return(ifelse(test = any(q$f!=1), yes = FALSE, no = TRUE))
    }
    else {
      q <- u[duplicated(u[,c("individual_id", "condition")])==FALSE,]
      q$f <- 1
      q <- aggregate(f~individual_id+condition, 
                     data = q, FUN = sum, drop = FALSE)
      q$f[is.null(q$f)|is.na(q$f)] <- 0
      return(ifelse(test = any(q$f!=1), yes = FALSE, no = TRUE))
    }
  }
  
  u$individual_org_name <- u$individual_id
  if("replicate" %in% colnames(u)) {
    
    has_replicates <- TRUE
    
    u$replicate_org_name <- u$replicate
    u$replicate_id <- u$replicate
    u$sample_id <- paste0(u$individual_id, '|', u$condition, '|', u$replicate)
    u$sample_id <- as.numeric(as.factor(u$sample_id))
    
    m <- u[duplicated(u$sample_id)==FALSE, c("sample_id",
                                             "individual_id", 
                                             "individual_org_name",
                                             "condition", 
                                             "replicate_org_name",
                                             "replicate_id")]
    
    # if replicate column is provided BUT only one replicate is available
    # per individual -> do analysis without replicates
    k <- u[duplicated(u[,c("individual_id","condition","replicate_id")])==FALSE,
           c("individual_id")]
    if(all(table(k)==1)) {
      has_replicates <- FALSE
    }
  } else {
    
    has_replicates <- FALSE
    
    u$sample_id <- paste0(u$individual_id, '|', u$condition)
    u$sample_id <- as.numeric(as.factor(u$sample_id))
    # u$individual_id <- paste0(u$individual_id, '|', u$condition)
    
    m <- u[duplicated(u$sample_id)==FALSE, c("sample_id", 
                                             "individual_id",
                                             "individual_org_name",
                                             "condition")]
  }
  
  u <- complete(u[,c("sample_id", "gene_name", "gene_usage_count")], 
                sample_id, gene_name, fill = list(gene_usage_count=0))
  
  n <- aggregate(gene_usage_count~sample_id, data = u, FUN = sum)
  n$total_usage_count <- n$gene_usage_count
  n$gene_usage_count <- NULL
  
  u <- merge(x = u, y = m, by = "sample_id")
  u <- merge(x = u, y = n, by = "sample_id")
  u$gene_usage_prop <- u$gene_usage_count/u$total_usage_count
  
  # get Y matrix
  Y <- acast(data = u, 
             formula = gene_name~sample_id,
             drop = FALSE, 
             value.var = "gene_usage_count",
             fill = 0, 
             fun.aggregate = sum)
  sample_ids <- colnames(Y)
  gene_names <- rownames(Y)
  
  N <- apply(X = Y, MARGIN = 2, FUN = sum)
  N <- as.numeric(N)
  
  # individual data
  individual_names <- character(length = length(sample_ids))
  individual_org_names <- character(length = length(sample_ids))
  
  # replicates data
  replicate_names <- character(length = length(sample_ids))
  replicate_org_names <- character(length = length(sample_ids))
  
  for(i in 1:length(sample_ids)) {
    individual_names[i] <- m$individual_id[m$sample_id == sample_ids[i]][1]
    individual_org_names[i] <- m$individual_org_name[
      m$sample_id == sample_ids[i]][1]
    
    if(has_replicates) {
      replicate_names[i] <- m$replicate_id[m$sample_id == sample_ids[i]][1]
      replicate_org_names[i] <- m$replicate_org_name[
        m$sample_id == sample_ids[i]][1]
    }
  }
  individual_ids <- as.numeric(as.factor(individual_names))
  replicate_ids <- as.numeric(as.factor(replicate_names))
  tr <- table(paste0(individual_ids, '|', replicate_ids))
  has_balanced_replicates <- ifelse(test = all(tr==tr[1]), yes=TRUE, no=FALSE)
  
  cos <- get_condition_of_sample(sample_ids = sample_ids, m = m) 
  coi <- get_condition_of_individual(individual_ids = individual_ids,
                                     individual_names = individual_names,
                                     m = m)
  
  has_conditions <- max(cos$condition_ids)>1
  
  has_paired_data <- check_paired(
    u = u, 
    has_balanced_replicates = has_balanced_replicates, 
    has_replicates = has_replicates, 
    has_conditions = has_conditions)
  
  return(list(Y = Y, 
              N = as.array(as.numeric(N)), 
              N_sample = ncol(Y), 
              N_gene = nrow(Y),
              gene_names = gene_names,
              sample_names = sample_ids, 
              condition_id_of_sample = cos$condition_ids,
              condition_name_of_sample = cos$condition_names,
              condition_id_of_individual = coi$condition_ids,
              condition_name_of_individual = coi$condition_names,
              N_condition = max(cos$condition_ids),
              individual_id = individual_ids,
              individual_names = individual_names,
              individual_org_names = individual_org_names,
              N_individual = max(individual_ids),
              replicate_id = replicate_ids,
              replicate_names = replicate_names,
              replicate_org_names = replicate_org_names,
              N_replicate = max(replicate_ids),
              proc_ud = u,
              has_replicates = has_replicates,
              has_conditions = has_conditions,
              has_balanced_replicates = has_balanced_replicates,
              has_paired_data = has_paired_data))
}

# Description:
# get the appropriate model
get_model <- function(has_replicates,
                      has_conditions,
                      has_balanced_replicates,
                      has_paired_data,
                      paired) {
  
  model_type <- ifelse(test = has_conditions, yes = "DGU", no = "GU")
  
  if(paired == TRUE & has_balanced_replicates == FALSE) {
    stop("For paired analysis with replicates, you need balanced replicates")
  }
  if(paired == TRUE & has_paired_data == FALSE) {
    stop("paired data selected, but data is not paired")
  }
  
  if(model_type == "GU") {
    if(has_replicates) {
      model <- stanmodels$gu_rep
      pars <- c("phi", "kappa", 
                "sigma_individual", "sigma_beta_rep",
                "beta_sample", "beta_individual", "beta_condition",
                "Yhat_rep", "Yhat_rep_prop", "Yhat_condition_prop", 
                "log_lik", "theta")
      model_name <- "GU_rep"
    } 
    else {
      model <- stanmodels$gu
      pars <- c("phi", "kappa", 
                "sigma_individual",
                "beta_individual", "beta_condition",
                "Yhat_rep", "Yhat_rep_prop", "Yhat_condition_prop", 
                "log_lik", "theta")
      model_name <- "GU"
    }
  } else {
    if(has_replicates) {
      model <- stanmodels$dgu_rep
      pars <- c("phi", "kappa", "alpha", 
                "sigma_condition", "sigma_individual", "sigma_beta_rep",
                "beta_sample", "beta_individual", "beta_condition", 
                "Yhat_rep", "Yhat_rep_prop", "Yhat_condition_prop", 
                "log_lik", "dgu", "dgu_prob", "theta")
      model_name <- "DGU_rep"
      if(paired) {
        model <- stanmodels$dgu_paired_rep
        pars <- c("phi", "kappa", "alpha", 
                  "sigma_condition", "sigma_individual", "sigma_beta_rep",
                  "mu_individual", "beta_sample", "beta_individual", 
                  "beta_condition", 
                  "Yhat_rep", "Yhat_rep_prop", "Yhat_condition_prop", 
                  "log_lik", "dgu", "dgu_prob", "theta")
        model_name <- "DGU_paired_rep"
      }
    } 
    else {
      model <- stanmodels$dgu
      pars <- c("phi", "kappa", "alpha", 
                "sigma_condition", "sigma_individual",
                "beta_individual", "beta_condition", 
                "Yhat_rep", "Yhat_rep_prop", "Yhat_condition_prop", 
                "log_lik", "dgu", "dgu_prob", "theta")
      model_name <- "DGU"
      if(paired) {
        model <- stanmodels$dgu_paired
        pars <- c("phi", "kappa", "alpha", 
                  "sigma_condition", "sigma_individual",
                  "mu_individual", "beta_individual", "beta_condition",
                  "Yhat_rep", "Yhat_rep_prop", "Yhat_condition_prop", 
                  "log_lik", "dgu", "dgu_prob", "theta")
        model_name <- "DGU_paired"
      }
    }
  }
  
  return(list(model = model, 
              model_name = model_name,
              model_type = model_type,
              pars = pars,
              has_replicates = has_replicates,
              has_conditions = has_conditions,
              has_balanced_replicates = has_balanced_replicates,
              has_paired_data = has_paired_data,
              paired = paired))
}

# Description:
# get the appropriate model
get_model_debug <- function(has_replicates, 
                            has_conditions, 
                            has_balanced_replicates) {
  model_type <- ifelse(test = has_conditions, yes = "DGU", no = "GU")
  
  if(model_type == "GU") {
    if(has_replicates) {
      model <- rstan::stan_model(file = "inst/stan/gu_rep.stan")
      pars <- c("phi", "kappa", 
                "sigma_individual", "sigma_beta_rep",
                "beta_sample", "beta_individual", "beta_condition",
                "Yhat_rep", "Yhat_rep_prop", "Yhat_condition_prop", 
                "log_lik", "theta")
      model_name <- "GU_rep"
    } 
    else {
      model <- rstan::stan_model(file = "inst/stan/gu.stan")
      pars <- c("phi", "kappa", 
                "sigma_individual",
                "beta_individual", "beta_condition",
                "Yhat_rep", "Yhat_rep_prop", "Yhat_condition_prop", 
                "log_lik", "theta")
      model_name <- "GU"
    }
  } 
  else {
    if(has_replicates) {
      model <- rstan::stan_model(file = "inst/stan/dgu_rep.stan")
      pars <- c("phi", "kappa", "alpha", 
                "sigma_condition", "sigma_individual", "sigma_beta_rep",
                "beta_sample", "beta_individual", "beta_condition", 
                "Yhat_rep", "Yhat_rep_prop", "Yhat_condition_prop", 
                "log_lik", "dgu", "dgu_prob", "theta")
      model_name <- "DGU_rep"
    } 
    else {
      model <- rstan::stan_model(file = "inst/stan/dgu.stan")
      pars <- c("phi", "kappa", "alpha", 
                "sigma_condition", "sigma_individual",
                "beta_individual", "beta_condition", 
                "Yhat_rep", "Yhat_rep_prop", "Yhat_condition_prop", 
                "log_lik", "dgu", "dgu_prob", "theta")
      model_name <- "DGU"
    }
  }
  
  return(list(model = model, 
              model_name = model_name,
              model_type = model_type,
              pars = pars,
              has_replicates = has_replicates,
              has_conditions = has_conditions))
}
