
# Description
# Parse input data for GU analysis
get_gu_usage <- function(u) {
  k <- base::paste0(u$sample_id, '|', u$condition)
  u$sample_id <- k
  
  # makes sure that we add usage = 0 for missing genes of certain samples 
  u <- tidyr::complete(u, sample_id, gene_name,
                       fill = base::list(gene_usage_count = 0))
  u$condition <- base::do.call(rbind, base::strsplit(
    x = u$sample_id, split = "\\|"))[,2]
  
  # format data
  n <- stats::aggregate(gene_usage_count~sample_id+condition, 
                        FUN = base::sum, data = u)
  n$total_usage_count <- n$gene_usage_count
  n$gene_usage_count <- NULL
  
  u <- base::merge(x = u, y = n, by = c("sample_id", "condition"), all.x = TRUE)
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

