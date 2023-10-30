
# Description
# Parse input data for GU analysis
get_gu_usage <- function(u) {
  k <- paste0(u$sample_id, '|', u$condition)
  u$sample_id <- k
  
  # makes sure that we add usage = 0 for missing genes of certain samples 
  u <- tidyr::complete(u, sample_id, gene_name,
                       fill = list(gene_usage_count = 0))
  u$condition <- do.call(rbind, strsplit(
    x = u$sample_id, split = "\\|"))[,2]
  
  # format data
  n <- aggregate(gene_usage_count~sample_id+condition, FUN = sum, data = u)
  n$total_usage_count <- n$gene_usage_count
  n$gene_usage_count <- NULL
  
  u <- merge(x = u, y = n, by = c("sample_id", "condition"), all.x = TRUE)
  u$gene_usage_prop <- u$gene_usage_count/u$total_usage_count
  
  # get Y matrix
  Y <- reshape2::acast(data = u, 
                       formula = gene_name~sample_id,
                       drop = FALSE, 
                       value.var = "gene_usage_count",
                       fill = 0, 
                       fun.aggregate = sum)
  sample_ids <- colnames(Y)
  gene_names <- rownames(Y)
  
  N <- apply(X = Y, MARGIN = 2, FUN = sum)
  N <- as.numeric(N)
  
  group_names <- character(length = length(sample_ids))
  for(i in seq_len(length(sample_ids))) {
    group_names[i] <- u$condition[u$sample == sample_ids[i]][1]
  }
  group_id <- as.numeric(as.factor(group_names))
  
  return (list(Y = Y, 
               N = as.numeric(N), 
               N_sample = ncol(Y), 
               N_gene = nrow(Y),
               gene_names = gene_names,
               sample_names = sample_ids, 
               group_id = group_id,
               group_names = group_names,
               N_group = max(group_id),
               proc_ud = u))
}

