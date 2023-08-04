# Description:
# From the input data format guess the type of the analysis: paired or unpaired
get_analysis_type <- function(u) {
  if(base::all(c("gene_usage_count_1","gene_usage_count_2")%in%
               base::colnames(u))){
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
# Parse input data for unpaired DGU analysis
get_unpaired_usage <- function(u) {
  k <- base::paste0(u$sample_id, '_', u$condition)
  u$sample_id <- k
  
  # makes sure that we add usage = 0 for missing genes of certain samples 
  u <- tidyr::complete(u, sample_id, gene_name, 
                       fill = base::list(gene_usage_count = 0))
  u$condition<-do.call(rbind, base::strsplit(x = u$sample_id, split="\\_"))[,2]
  
  # format data
  n <- stats::aggregate(gene_usage_count~sample_id+condition, 
                        FUN = base::sum, data = u)
  n$total_usage_count <- n$gene_usage_count
  n$gene_usage_count <- NULL
  
  u <- base::merge(x = u, y = n, by = c("sample_id", "condition"))
  u$gene_usage_prop <- u$gene_usage_count/u$total_usage_count
  
  # get Y matrix
  Y <- reshape2::acast(data = u, 
                       formula = gene_name~sample_id,
                       drop = FALSE, 
                       value.var = "gene_usage_count",
                       fill = 0, 
                       fun.aggregate = base::sum)
  sample_ids <- base::colnames(Y)
  gene_names <- base::rownames(Y)
  
  N <- base::apply(X = Y, MARGIN = 2, FUN = base::sum)
  N <- base::as.numeric(N)
  
  X <- base::numeric(length = base::length(N))
  
  # sorted unique conditions
  cs <- base::sort(base::unique(u$condition), decreasing = TRUE)
  # design vector X
  X <- base::ifelse(test = sample_ids %in% u$sample_id[u$condition == cs[1]], 
                    yes = +1, no = -1)
  
  # group
  group_names <- base::character(length = base::length(sample_ids))
  for(i in base::seq_len(base::length(sample_ids))) {
    group_names[i] <- u$condition[u$sample == sample_ids[i]][1]
  }
  group_id <- base::as.numeric(base::as.factor(group_names))
  
  # add contrast
  contrast <- base::paste0(base::unique(group_names[group_id == 1][1]), " - ", 
                           base::unique(group_names[group_id == 2][1]))
  
  return (base::list(Y = Y, 
                     N = base::as.numeric(N), 
                     N_sample = base::ncol(Y), 
                     N_gene = base::nrow(Y),
                     X = X,
                     gene_names = gene_names,
                     sample_names = sample_ids,
                     group_id = group_id,
                     group_names = group_names,
                     N_group = base::max(group_id),
                     proc_ud = u,
                     contrast = contrast))
}


# Description:
# Parse input data for paired DGU analysis
get_paired_usage <- function(u) {
  
  # makes sure that we add usage = 0 for missing genes of certain samples 
  u <- tidyr::complete(u, sample_id, gene_name, 
                       fill = base::list(gene_usage_count_1 = 0,
                                         gene_usage_count_2 = 0))
  
  # format data
  u_1 <- stats::aggregate(gene_usage_count_1~sample_id, 
                         FUN = base::sum, data = u)
  u_1$total_usage_count_1 <- u_1$gene_usage_count_1
  u_1$gene_usage_count_1 <- NULL
  
  u_2 <- stats::aggregate(gene_usage_count_2~sample_id, 
                         FUN = base::sum, data = u)
  u_2$total_usage_count_2 <- u_2$gene_usage_count_2
  u_2$gene_usage_count_2 <- NULL
  
  u <- base::merge(x = u, y = u_1, by = c("sample_id"))
  u <- base::merge(x = u, y = u_2, by = c("sample_id"))
  u$gene_usage_prop_1 <- u$gene_usage_count_1/u$total_usage_count_1
  u$gene_usage_prop_2 <- u$gene_usage_count_2/u$total_usage_count_2
  
  # get matrices
  # Y_1
  Y_1 <- reshape2::acast(data = u, 
                         formula = gene_name~sample_id,
                         drop = FALSE, 
                         value.var = "gene_usage_count_1",
                         fill = 0, 
                         fun.aggregate = base::sum)
  Y_1 <- Y_1[base::rownames(Y_1),]
  
  
  # Y_2
  Y_2 <- reshape2::acast(data = u, 
                         formula = gene_name~sample_id,
                         drop = FALSE, 
                         value.var = "gene_usage_count_2",
                         fill = 0, 
                         fun.aggregate = base::sum)
  # order Y_1 and Y_2 by same gene names
  Y_2 <- Y_2[base::rownames(Y_1),] 
  
  
  # compute total usage
  N_1 <- base::apply(X = Y_1, MARGIN = 2, FUN = base::sum)
  N_2 <- base::apply(X = Y_2, MARGIN = 2, FUN = base::sum)
  
  
  return (base::list(Y_1 = Y_1, 
                     Y_2 = Y_2, 
                     N_1 = N_1,
                     N_2 = N_2,
                     N_sample = base::ncol(Y_1), 
                     N_gene = base::nrow(Y_1), 
                     N_group = 1,
                     gene_names = base::rownames(Y_1),
                     sample_names = base::colnames(Y_1), 
                     proc_ud = u))
}


# Description
# Parse input data for GU analysis
get_gu_usage <- function(u) {
  k <- base::paste0(u$sample_id, '_', u$condition)
  u$sample_id <- k
  
  # makes sure that we add usage = 0 for missing genes of certain samples 
  u <- tidyr::complete(u, sample_id, gene_name,
                       fill = base::list(gene_usage_count = 0))
  u$condition <- base::do.call(rbind, base::strsplit(
    x = u$sample_id, split = "\\_"))[,2]
  
  # format data
  n <- stats::aggregate(gene_usage_count~sample_id+condition, 
                        FUN = base::sum, data = u)
  n$total_usage_count <- n$gene_usage_count
  n$gene_usage_count <- NULL
  
  u <- base::merge(x = u, y = n, by = c("sample_id", "condition"), all.x = T)
  u$gene_usage_prop <- u$gene_usage_count/u$total_usage_count
  
  # get Y matrix
  Y <- reshape2::acast(data = u, 
                       formula = gene_name~sample_id,
                       drop = FALSE, 
                       value.var = "gene_usage_count",
                       fill = 0, 
                       fun.aggregate = base::sum)
  sample_ids <- base::colnames(Y)
  gene_names <- base::rownames(Y)
  
  N <- base::apply(X = Y, MARGIN = 2, FUN = base::sum)
  N <- base::as.numeric(N)
  
  group_names <- base::character(length = base::length(sample_ids))
  for(i in base::seq_len(base::length(sample_ids))) {
    group_names[i] <- u$condition[u$sample == sample_ids[i]][1]
  }
  group_id <- base::as.numeric(base::as.factor(group_names))
  
  return (base::list(Y = Y, 
                     N = base::as.numeric(N), 
                     N_sample = base::ncol(Y), 
                     N_gene = base::nrow(Y),
                     gene_names = gene_names,
                     sample_names = sample_ids, 
                     group_id = group_id,
                     group_names = group_names,
                     N_group = base::max(group_id),
                     proc_ud = u))
}

