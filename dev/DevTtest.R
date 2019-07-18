data("IGHV_HCV", package = "IgGeneUsage")
N <- aggregate(gene_usage_count~sample_id, FUN = sum, data = IGHV_HCV)
N$total_usage_count <- N$gene_usage_count
N$gene_usage_count <- NULL

source("R/Usage.R")
source("R/Util.R")

rm(viz)



# two sided t.test
getDebugTtest <- function(usage.data) {
  getTTest <- function(x, Ys, Xs, Ns) {
    return(try(stats::t.test(Ys[x, ]/Ns*100~Xs)))
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

  tout <- suppressWarnings(expr = lapply(X = 1:usage.data$N_gene,
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


samples <- unique(IGHV_HCV$sample_id)

complete.t.test.stats <- c()
for(s in 1:length(samples)) {
  usage.data <- getUsageData(usage = IGHV_HCV[!IGHV_HCV$sample_id %in% samples[s], ])
  t.test.stats <- getDebugTtest(usage.data = usage.data)
  t.test.stats$boot.id <- s
  complete.t.test.stats <- rbind(complete.t.test.stats, t.test.stats)
  rm(t.test.stats, usage.data)
}

boxplot(complete.t.test.stats$t.test.bonf.pvalue[complete.t.test.stats$gene_name == "IGHV1-58"])
x <- complete.t.test.stats[complete.t.test.stats$gene_name == "IGHV1-58", ]
y <- IGHV_HCV[IGHV_HCV$gene_name == "IGHV1-58", ]
y <- merge(x = y, y = N, by = c("sample_id"))
y$gene_usage_p <- y$gene_usage_count/y$total_usage_count

plot(y$gene_usage_p)


# frequentist tests, merge data
t.test.stats <- getTTestStats(usage.data = usage.data)
u.test.stats <- getManUStats(usage.data = usage.data)
test.stats <- merge(x = t.test.stats, y = u.test.stats, by = "gene_name")
