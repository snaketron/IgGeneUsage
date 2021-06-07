


# Description:
# Parse input data
getUsageData <- function(usage) {
  # if the same sample_id is present in both conditions
  key <- paste(usage$sample_id, usage$condition, sep = '_')
  if(length(unique(key)) != length(unique(usage$sample_id))) {
    warning("Same sample_id in both conditions, sample_id's extra coded")
    usage$sample_id <- key
    rm(key)
  }

  # get Y data, fill empty combinations with 0
  Y <- reshape2::acast(data = usage, formula = gene_name~sample_id,
                       drop = FALSE, value.var = "gene_usage_count",
                       fill = 0, fun.aggregate = sum)

  sample_ids <- colnames(Y)
  gene_names <- rownames(Y)

  # get N data
  N.data <- stats::aggregate(gene_usage_count~sample_id, data = usage,
                             FUN = sum, drop = FALSE)
  N <- N.data$gene_usage_count
  names(N) <- N.data$sample_id
  rm(N.data)
  N <- N[sample_ids]

  # get X data
  usage <-  usage[usage$sample_id %in% sample_ids, ]
  usage <- usage[duplicated(usage[, c("sample_id")]) == FALSE, ]

  X.data <- stats::aggregate(condition~sample_id, data = usage, FUN = unique)
  X <- X.data$condition
  names(X) <- X.data$sample_id
  rm(X.data)
  X <- X[sample_ids]

  # get X mapping
  X.unique <- sort(x = unique(X), decreasing = TRUE)
  x1 <- which(X == X.unique[1])
  x2 <- which(X == X.unique[2])
  Xmap <- numeric(length = length(X))
  Xmap[x1] <- 1
  Xmap[x2] <- -1
  
  
  # process usage
  for(i in 1:ncol(Y)) {
    processed.usage.data <- data.frame(gene_usage_count = Y[,i],
                                       total_usage_count = sum(Y[,i]),
                                       gene_usage_prop = Y[,i]/sum(Y[,i]),
                                       sample_id = sample_ids[i],
                                       condition = X[i],
                                       stringsAsFactors = F)
    
  }

  return (list(Y = Y, N = N, N_sample = ncol(Y), N_gene = nrow(Y),
               X = Xmap, Xorg = X, gene_names = gene_names,
               sample_names = sample_ids,
               processed.usage.data = processed.usage.data))
}



# Description:
# Takes an input usage.data of class SummarizedExperiment,
# and creates data.frame. Given that this function is used
# before the input control, it also does some rudimentary checks
convertSummarizedExperiment <- function(usage.data.se) {
  # get count data
  count.data <- SummarizedExperiment::assay(x = usage.data.se)
  count.data <- reshape2::melt(
    data = count.data, value.name = "gene_usage_count", as.is = TRUE)
  colnames(count.data) <- c("gene_name", "sample_id", "gene_usage_count")


  # get the sample (repertoire) data and convert it to data.frame
  coldata <- SummarizedExperiment::colData(x = usage.data.se)
  coldata <- base::as.data.frame(coldata)
  if(nrow(coldata) == 0 | ncol(coldata) == 0) {
    stop("colData(usage.data) is empty")
  }

  if(all(colnames(coldata) %in% c("condition", "sample_id")) == FALSE) {
    stop("colData(usage.data) must have 2 columns: 'condition' and 'sample_id'")
  }

  # merge both datasets into the standard 4-column IgGeneUsage input
  usage.data <- base::merge(x = count.data,
                            y = coldata[, c("condition", "sample_id")],
                            by = "sample_id", all = TRUE)

  return (usage.data)
}






getPmax <- function(glm.ext) {

  getPmaxGene <- function(x, beta.data) {
    p <- beta.data[,x]
    l <- length(p)
    o <- max(sum(p < 0)/l, sum(p>0)/l)
    return(o)
  }

  beta.data <- glm.ext$beta_gene
  pmax <- vapply(X = seq_len(length.out = ncol(beta.data)),
                 FUN = getPmaxGene,
                 beta.data = beta.data,
                 FUN.VALUE = numeric(length = 1))

  # convert from [0.5, 1] -> [0, 1]
  pmax <- 2 * pmax - 1
  return(pmax)
}




# Description
# Posterior-predictive check within repertoires
getPpcRepertoire <- function(glm,
                             usage.data,
                             hdi.level) {

  # summaries
  yhat <- summary(object = glm, pars = c("Yhat", "Yhat_individual"),
                  prob = c(0.5, (1-hdi.level)/2, 1-(1-hdi.level)/2))
  yhat <- yhat$summary
  yhat <- data.frame(yhat)
  colnames(yhat) <- c("mean", "se_mean", "sd", "median",
                      "L", "H", "Neff", "Rhat")
  yhat[, c("Rhat", "Neff")] <- NULL
  yhat$par <- rownames(yhat)
  yhat$par.name <- do.call(rbind, strsplit(x = yhat$par, split = "\\["))[, 1]
  par.indices <- do.call(rbind, strsplit(x = yhat$par, split = "\\["))[, 2]
  par.indices <- gsub(pattern = "\\]", replacement = '', x = par.indices)
  par.indices <- do.call(rbind, strsplit(x = par.indices, split = ','))
  class(par.indices) <- "numeric"
  yhat$G <- par.indices[, 1]
  yhat$R <- par.indices[, 2]
  rm(par.indices)

  # split
  yhat.count <- yhat[yhat$par.name == "Yhat", ]
  yhat.prop <- yhat[yhat$par.name == "Yhat_individual", ]
  rm(yhat)

  yhat.count$par <- NULL
  yhat.count$par.name <- NULL
  colnames(yhat.count)[1:6] <- paste("ppc", colnames(yhat.count)[1:6],
                                     "count", sep = "_")

  yhat.prop$par <- NULL
  yhat.prop$par.name <- NULL
  colnames(yhat.prop)[1:6] <- paste("ppc", colnames(yhat.prop)[1:6],
                                   "prop", sep = "_")

  yhat <- merge(x = yhat.count, y = yhat.prop, by = c("G", "R"))
  rm(yhat.count, yhat.prop)

  yhat$condition <- NA
  yhat$sample_name <- NA
  yhat$gene_name <- NA
  yhat$observed_count <- NA
  yhat$observed_prop <- NA

  for(i in 1:nrow(yhat)) {
    yhat$sample_name[i] <- usage.data$sample_names[yhat$R[i]]
    yhat$gene_name[i] <- usage.data$gene_names[yhat$G[i]]
    yhat$observed_count[i] <- usage.data$Y[yhat$G[i], yhat$R[i]]
    yhat$observed_prop[i] <- yhat$observed_count[i]/usage.data$N[yhat$R[i]]
    yhat$condition[i] <- usage.data$Xorg[yhat$sample_name[i]]
  }

  return (yhat)
}



# Description
# Posterior-predictive check across genes
getPpcGene <- function(glm,
                       usage.data,
                       hdi.level) {

  # summaries
  yhat <- summary(object = glm, pars = c("Yhat_gene"),
                  prob = c(0.5, (1-hdi.level)/2, 1-(1-hdi.level)/2))
  yhat <- yhat$summary
  yhat <- data.frame(yhat)
  colnames(yhat) <- c("mean", "se_mean", "sd", "median",
                      "L", "H", "Neff", "Rhat")
  yhat[, c("Rhat", "Neff")] <- NULL
  yhat$par <- rownames(yhat)
  yhat$par.name <- do.call(rbind, strsplit(x = yhat$par, split = "\\["))[, 1]
  par.indices <- do.call(rbind, strsplit(x = yhat$par, split = "\\["))[, 2]
  par.indices <- gsub(pattern = "\\]", replacement = '', x = par.indices)
  par.indices <- do.call(rbind, strsplit(x = par.indices, split = ','))
  class(par.indices) <- "numeric"
  yhat$X <- par.indices[, 1]
  yhat$G <- par.indices[, 2]
  rm(par.indices)

  yhat$par <- NULL
  yhat$par.name <- NULL
  colnames(yhat)[1:6] <- paste("ppc", colnames(yhat)[1:6], "prop", sep = "_")

  yhat$X <- ifelse(test = yhat$X == 1, yes = 1, no = -1)
  yhat$condition <- NA
  yhat$gene_name <- NA
  yhat$observed_prop <- NA
  for(i in 1:nrow(yhat)) {
    yhat$gene_name[i] <- usage.data$gene_names[yhat$G[i]]
    # # TODO: here check
    # browser()
    yhat$observed_prop[i] <- mean(usage.data$Y[yhat$G[i], ]/usage.data$N)
    yhat$condition[i] <- usage.data$Xorg[usage.data$X == yhat$X[i]][1]
  }

  return (yhat)
}




# two sided t.test
getTTestStats <- function(usage.data) {
  getTTest <- function(x, Ys, Xs, Ns) {
    return(try(stats::t.test((Ys[x, ]/Ns)~Xs)))
  }
  getTTestSummary <- function(x) {
    if(inherits(x = x, what = 'try-error') == TRUE) {
      return(data.frame(t.test.pvalue = NA,
                        t.test.tvalue = NA,
                        t.test.L95 = NA,
                        t.test.H95 = NA,
                        stringsAsFactors = FALSE))
    }
    return(data.frame(t.test.pvalue = x$p.value,
                      t.test.tvalue = x$statistic,
                      t.test.L95 = x$conf.int[1],
                      t.test.H95 = x$conf.int[2],
                      stringsAsFactors = FALSE))
  }

  tout <- suppressWarnings(
    expr = lapply(X = seq_len(length.out = usage.data$N_gene),
                  FUN = getTTest,
                  Ys = usage.data$Y,
                  Xs = usage.data$X,
                  Ns = usage.data$N))

  tout.summary <- do.call(rbind, lapply(tout, getTTestSummary))
  tout.summary$gene_name <- usage.data$gene_names

  # multiple correction
  tout.summary$t.test.fdr.pvalue <- stats::p.adjust(
    p = tout.summary$t.test.pvalue, method = "fdr")

  return (tout.summary)
}





getManUStats <- function(usage.data) {

  getMTest <- function(x, Ys, Xs, Ns) {
    return(try(stats::wilcox.test((Ys[x, ]/Ns)~Xs)))
  }

  getMSummary <- function(x) {
    if(inherits(x = x, what = "try-error") == TRUE) {
      return(data.frame(u.test.pvalue = NA,
                        u.test.wvalue = NA,
                        stringsAsFactors = FALSE))
    }
    return(data.frame(u.test.pvalue = x$p.value,
                      u.test.wvalue = x$statistic,
                      stringsAsFactors = FALSE))
  }

  mout <- suppressWarnings(
    expr = lapply(X = seq_len(length.out = usage.data$N_gene),
                  FUN = getMTest,
                  Ys = usage.data$Y,
                  Xs = usage.data$X,
                  Ns = usage.data$N))

  mout.summary <- do.call(rbind, lapply(mout, getMSummary))
  mout.summary$gene_name <- usage.data$gene_names

  # multiple correction
  mout.summary$u.test.fdr.pvalue <- stats::p.adjust(
    p = mout.summary$u.test.pvalue, method = "fdr")

  return (mout.summary)
}


