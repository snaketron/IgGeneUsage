


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
  
  return (list(Y = Y, N = N, N_sample = ncol(Y), N_gene = nrow(Y),
               X = Xmap, Xorg = X, gene_names = gene_names,
               sample_names = sample_ids))
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
  usage.data <- base::merge(x = count.data, y = coldata, 
                            by = "sample_id", all = TRUE)
  
  return (usage.data)
}




# Description:
# Computes HDI given a vector, taken "Doing Bayesian Analysis"
getHdi <- function(vec, hdi.level) {
  sortedPts <- sort(vec)
  ciIdxInc <- floor(hdi.level * length(sortedPts))
  nCIs = length(sortedPts) - ciIdxInc
  ciWidth = rep(0 , nCIs)
  for (i in seq_len(length.out = nCIs)) {
    ciWidth[i] = sortedPts[i + ciIdxInc] - sortedPts[i]
  }
  HDImin = sortedPts[which.min(ciWidth)]
  HDImax = sortedPts[which.min(ciWidth) + ciIdxInc]
  HDIlim = c(HDImin, HDImax)
  return(HDIlim)
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

  return(pmax)
}





getPpc <- function(glm.ext,
                   usage.data,
                   hdi.level) {

  yhat.pct <- glm.ext$Yhat_individual
  yhat.count <- glm.ext$Yhat

  ppc <- c()
  for(i in seq_len(length.out = usage.data$N_sample)) {
    for(j in seq_len(length.out = usage.data$N_gene)) {
      yhat.i.count <- yhat.count[,j,i]
      hdi.count <- getHdi(vec = yhat.i.count, hdi.level = hdi.level)

      yhat.i.pct <- yhat.pct[,j,i]
      hdi.pct <- getHdi(vec = yhat.i.pct, hdi.level = hdi.level)

      # errors
      error.count <- abs(usage.data$Y[j,i]-yhat.i.count)
      error.pct <- abs(usage.data$Y[j,i]/usage.data$N[i]*100-yhat.i.pct)

      row <- data.frame(sample_id = usage.data$sample_names[i],
                        gene_name = usage.data$gene_names[j],
                        condition = usage.data$Xorg[i],
                        observed.count = usage.data$Y[j,i],
                        observed.pct = usage.data$Y[j,i]/usage.data$N[i]*100,
                        ppc.count.mean = mean(yhat.i.count),
                        ppc.count.median = stats::median(yhat.i.count),
                        ppc.count.L = hdi.count[1],
                        ppc.count.H = hdi.count[2],
                        ppc.pct.mean = mean(yhat.i.pct),
                        ppc.pct.median = stats::median(yhat.i.pct),
                        ppc.pct.L = hdi.pct[1],
                        ppc.pct.H = hdi.pct[2],
                        error.count.mean = mean(error.count),
                        error.count.median = stats::median(error.count),
                        error.pct.mean = mean(error.pct),
                        error.pct.median = stats::median(error.pct),
                        stringsAsFactors = FALSE)

      ppc <- rbind(ppc, row)
    }
  }

  return (ppc)
}




getGroupStats <- function(glm.ext, usage.data, hdi.level) {
  
  getGroupYhat <- function(x, gene.names, conditions, yhat.gene, 
                           hdi.level, usage.data) {

    hdi.1 <- getHdi(vec = yhat.gene[,1,x], hdi.level = hdi.level)
    hdi.2 <- getHdi(vec = yhat.gene[,2,x], hdi.level = hdi.level)

    # get mean count data
    x.1 <- which(usage.data$Xorg == conditions[1])
    real.pct.1 <- 0
    if(length(x.1) != 0) {
      real.pct.1 <- usage.data$Y[gene.names[x], x.1]/usage.data$N[x.1]*100
    }
    x.2 <- which(usage.data$Xorg == conditions[2])
    real.pct.2 <- 0
    if(length(x.2) != 0) {
      real.pct.2 <- usage.data$Y[gene.names[x], x.2]/usage.data$N[x.2]*100
    }

    # errors
    error.pct.1 <- abs(mean(real.pct.1)-yhat.gene[,1,x])
    error.pct.2 <- abs(mean(real.pct.2)-yhat.gene[,2,x])
    return(rbind(data.frame(gene_name = gene.names[x],
                            observed.mean = mean(real.pct.1),
                            ppc.mean = mean(yhat.gene[,1,x]),
                            ppc.L = hdi.1[1], ppc.H = hdi.1[2],
                            error.mean = mean(error.pct.1),
                            condition = conditions[1],
                            stringsAsFactors = FALSE),
                 data.frame(gene_name = gene.names[x],
                            observed.mean = mean(real.pct.2),
                            ppc.mean = mean(yhat.gene[,2,x]),
                            ppc.L = hdi.2[1], ppc.H = hdi.2[2],
                            error.mean = mean(error.pct.2),
                            condition = conditions[2],
                            stringsAsFactors = FALSE)))
  }

  conditions = c(unique(usage.data$Xorg[usage.data$X == 1]),
                 unique(usage.data$Xorg[usage.data$X == -1]))
  group.ppc <- lapply(X = seq_len(length.out = usage.data$N_gene),
                      FUN = getGroupYhat,
                      gene.names = usage.data$gene_names,
                      condition = conditions,
                      yhat.gene = glm.ext$Yhat_gene,
                      hdi.level = hdi.level,
                      usage.data = usage.data)
  return (do.call(rbind, group.ppc))
}




# two sided t.test
getTTestStats <- function(usage.data) {
  getTTest <- function(x, Ys, Xs, Ns) {
    return(try(stats::t.test((Ys[x, ]/Ns*100)~Xs)))
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
  tout.summary$t.test.bonf.pvalue <- stats::p.adjust(
    p = tout.summary$t.test.pvalue, method = "bonferroni")



  return (tout.summary)
}





getManUStats <- function(usage.data) {

  getMTest <- function(x, Ys, Xs, Ns) {
    return(try(stats::wilcox.test((Ys[x, ]/Ns*100)~Xs)))
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
  mout.summary$u.test.bonf.pvalue <- stats::p.adjust(
    p = mout.summary$u.test.pvalue, method = "bonferroni")

  return (mout.summary)
}


