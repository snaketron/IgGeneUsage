


# Description:
# Parse input data
getUsageData <- function(usage) {

  # get Y data
  # fill empty combinations with 0
  Y <- reshape2::acast(data = usage,
                       formula = gene_name~sample_id,
                       drop = FALSE,
                       value.var = "gene_usage_count",
                       fill = 0,
                       fun.aggregate = sum)
  sample_ids <- colnames(Y)
  gene_names <- rownames(Y)

  # get N data
  N.data <- aggregate(gene_usage_count~sample_id,
                      data = usage, FUN = sum,
                      drop = FALSE)
  N <- N.data$gene_usage_count
  names(N) <- N.data$sample_id
  rm(N.data)
  N <- N[sample_ids]


  # get X data
  usage <-  usage[usage$sample_id %in% sample_ids, ]
  usage <- usage[duplicated(usage[, c("sample_id")]) == FALSE, ]

  X.data <- aggregate(condition~sample_id,
                      data = usage, FUN = unique)
  X <- X.data$condition
  names(X) <- X.data$sample_id
  rm(X.data)
  X <- X[sample_ids]


  # get X mapping
  X.unique <- unique(X)
  x1 <- which(X == X.unique[1])
  x2 <- which(X == X.unique[2])
  Xmap <- numeric(length = length(X))
  Xmap[x1] <- 1
  Xmap[x2] <- -1
  if(runif(n = 1, min = 0, max = 1) >= 0.5) {
    Xmap[x1] <- -1
    Xmap[x2] <- 1
  }


  usage.data <- list(Y = Y,
                     N = N,
                     N_sample = ncol(Y),
                     N_gene = nrow(Y),
                     X = Xmap,
                     Xorg = X,
                     gene_names = gene_names,
                     sample_names = sample_ids)


  return (usage.data)
}






# Description:
# Check the optional (...) parameters, if provided.
checkDotParameters <- function(...) {

  checkAdaptDelta <- function(adapt_delta) {
    if(length(adapt_delta) != 1) {
      stop("adapt_delta must be in range (0, 1) (default = 0.8).")
    }

    if(is.numeric(adapt_delta) == FALSE) {
      stop("adapt_delta must be in range (0, 1)")
    }

    if(adapt_delta >= 1 | adapt_delta <= 0) {
      stop("adapt_delta must be in range (0, 1)")
    }
  }

  checkMaxTreedepth <- function(max_treedepth) {
    if(length(max_treedepth) != 1) {
      stop("max_treedepth is numeric parameter.")
    }

    if(is.numeric(max_treedepth) == FALSE) {
      stop("max_treedepth is numeric parameter.")
    }

    if(max_treedepth < 5) {
      stop("max_treedepth >= 5 (default = 10).")
    }
  }

  checkCvFold <- function(cv.fold) {
    if(length(cv.fold) != 1) {
      stop("cv.fold must be in range (0, 1) (default = 0.66).")
    }

    if(is.numeric(cv.fold) == FALSE) {
      stop("cv.fold must be in range (0, 1)")
    }

    if(cv.fold >= 1 | cv.fold <= 0) {
      stop("cv.fold must be in range (0, 1)")
    }
  }

  checkNtree <- function(ntree) {
    if(length(ntree) != 1) {
      stop("ntree is numeric parameter.")
    }

    if(is.numeric(ntree) == FALSE) {
      stop("ntree is numeric parameter.")
    }

    if(ntree < 100) {
      stop("ntree >= 100 (default = 500).")
    }
  }

  checkVerbose <- function(verbose) {
    if(length(verbose) != 1) {
      stop("verbose is a logical parameter.")
    }

    if(is.logical(verbose) == FALSE) {
      stop("verbose is a logical parameter.")
    }
  }

  checkRefresh <- function(refresh) {
    if(length(refresh) != 1) {
      stop("refresh is numeric parameter.")
    }

    if(is.numeric(refresh) == FALSE) {
      stop("refresh is a numeric parameter.")
    }

    return (refresh)
  }

  available.names <- c("adapt_delta",
                       "max_treedepth",
                       "ntree",
                       "cv.fold",
                       "refresh",
                       "verbose")
  default.values <- list(adapt_delta = 0.8,
                         max_treedepth = 10,
                         ntree = 1000,
                         cv.fold = 0.66,
                         refresh = 250,
                         verbose = TRUE)

  # get the optional parameters
  dot.names <- names(list(...))

  if(length(dot.names) > 0) {
    if(any(dot.names %in% available.names) == FALSE) {
      wrong.names <- dot.names[!dot.names %in% available.names]
      stop(paste("Unknown optional parameter were provided! The following
                 optional parameters are available:", dot.names, sep = ' '))
    }
  }

  # check each parameter
  for(p in dot.names) {
    if(is.null(list(...)[[p]]) || is.na(list(...)[[p]])) {
      stop(paste("optional parameter ", p, " can't be NULL", sep = ''))
    }
    if(p == "adapt_delta") {
      checkAdaptDelta(adapt_delta = list(...)[[p]])
      default.values[["adapt_delta"]] <- list(...)[[p]]
    }
    if(p == "max_treedepth") {
      checkMaxTreedepth(max_treedepth = list(...)[[p]])
      default.values[["max_treedepth"]] <- list(...)[[p]]
    }
    if(p == "cv.fold") {
      checkCvFold(cv.fold = list(...)[[p]])
      default.values[["cv.fold"]] <- list(...)[[p]]
    }
    if(p == "ntree") {
      checkNtree(ntree = list(...)[[p]])
      default.values[["ntree"]] <- list(...)[[p]]
    }
    if(p == "refresh") {
      checkRefresh(refresh = list(...)[[p]])
      default.values[["refresh"]] <- list(...)[[p]]
    }
    if(p == "verbose") {
      checkVerbose(verbose = list(...)[[p]])
      default.values[["verbose"]] <- list(...)[[p]]
    }
  }

  return (default.values)
}





# Description:
# Provided the input arguments, this function checks their validity. It
# stops the execution if a problem is encountered and prints out warnings.
checkInput <- function(usage.data,
                       mcmc.chains,
                       mcmc.cores,
                       mcmc.steps,
                       mcmc.warmup,
                       hdi.level) {

  checkUsageData <- function(usage.data) {
    if(class(usage.data) != "data.frame") {
      stop("usage.data must be data.frame")
    }

    if(nrow(usage.data) <= 1) {
      stop("usage.data must contain at least 2 rows (usage of one
           gene in two conditions.")
    }

    if(ncol(usage.data) != 4) {
      stop("usage.data must contain the following 4 columns: 'sample_id',
           'condition', 'gene_name' and 'gene_usage_count'")
    }

    correct.col <- c("sample_id", "condition", "gene_name", "gene_usage_count")
    if(all(colnames(usage.data) %in% correct.col) == FALSE) {
      stop("usage.data must contain the following columns: 'sample_id',
           'condition', 'gene_name' and 'gene_usage_count'")
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


  if(is.null(usage.data) | missing(usage.data) |
     is.null(mcmc.chains) | missing(mcmc.chains) |
     is.null(mcmc.steps) | missing(mcmc.steps) |
     is.null(mcmc.warmup) | missing(mcmc.warmup) |
     is.null(mcmc.cores) | missing(mcmc.cores) |
     is.null(hdi.level) | missing(hdi.level)) {
    stop("arguments must be non-NULL/specified")
  }


  checkUsageData(usage.data = usage.data)
  checkMcmcIterations(mcmc.steps = mcmc.steps)
  checkMcmcWarmup(mcmc.warmup = mcmc.warmup)
  checkMcmcChains(mcmc.chains = mcmc.chains)
  checkMcmcCores(mcmc.cores = mcmc.cores)
  checkHdi(hdi.level = hdi.level)
}





# Description:
# Computes HDI given a vector, taken "Doing Bayesian Analysis"
getHdi <- function(vec,
                   hdi.level) {
  sortedPts <- sort(vec)
  ciIdxInc <- floor(hdi.level * length(sortedPts))
  nCIs = length(sortedPts) - ciIdxInc
  ciWidth = rep(0 , nCIs)
  for (i in 1:nCIs) {
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
  pmax <- sapply(X = 1:ncol(beta.data),
                 FUN = getPmaxGene,
                 beta.data = beta.data)

  return(pmax)
}





getPpc <- function(glm.ext,
                   usage.data,
                   hdi.level) {
  yhat.pct <- glm.ext$Yhat_individual
  yhat.raw <- glm.ext$Yhat

  ppc <- c()
  for(i in 1:usage.data$N_sample) {
    for(j in 1:usage.data$N_gene) {
      yhat.i.raw <- yhat.raw[,j,i]
      hdi.raw <- getHdi(vec = yhat.i.raw, hdi.level = hdi.level)

      yhat.i.pct <- yhat.pct[,j,i]
      hdi.pct <- getHdi(vec = yhat.i.pct, hdi.level = hdi.level)


      row <- data.frame(sample_id = usage.data$sample_names[i],
                        gene_name = usage.data$gene_names[j],
                        condition = usage.data$Xorg[i],
                        raw = usage.data$Y[j,i],
                        raw.pct = usage.data$Y[j,i]/usage.data$N[i]*100,
                        ppc.raw.mean = mean(yhat.i.raw),
                        ppc.raw.median = median(yhat.i.raw),
                        ppc.raw.L = hdi.raw[1],
                        ppc.raw.H = hdi.raw[2],
                        ppc.pct.mean = mean(yhat.i.pct),
                        ppc.pct.median = median(yhat.i.pct),
                        ppc.pct.L = hdi.pct[1],
                        ppc.pct.H = hdi.pct[2])

      ppc <- rbind(ppc, row)
    }
  }

  return (ppc)
}





getBcStats <- function(glm.ext,
                       usage.data) {

  # Description:
  # Bhattacharyya Coefficient of two distribution
  # Taken from: source("http://tguillerme.github.io/R/bhatt.coef.R")
  getBCoef <- function(x, y, bw = bw.nrd0, ...) {
    #SANITIZING
    #x
    if(class(x) != 'numeric') {
      stop("'x' must be numeric.")
    }
    if(length(x) < 2) {
      stop("'x' need at least two data points.")
    }

    #y
    if(class(y) != 'numeric') {
      stop("'y' must be numeric.")
    }
    if(length(y) < 2) {
      stop("'y' need at least two data points.")
    }

    #bw
    if(length(bw) != 1) {
      stop("'bw' must be either a single numeric value or a single function.")
    }
    if(class(bw) != 'function') {
      if(class(bw) != 'numeric') {
        stop("'bw' must be either a single numeric value or a single function.")
      }
    }
    #Avoiding non-entire numbers
    if(class(bw) == 'numeric') {
      bw<-round(bw)
    }

    #BHATTACHARYYA COEFFICIENT
    #sum(sqrt(x relative counts in bin_i * y relative counts in bin_i))

    #Setting the right number of bins (i)
    if(class(bw) == 'function') {
      #Bin width
      band.width<-bw(c(x,y), ...)
      #Bin breaks
      bin.breaks<-seq(from=min(c(x,y)), to=max(c(x,y)+band.width), by=band.width) #adding an extra bandwith to the max to be sure to include all the data
      #Number of bins
      bin.n<-length(bin.breaks)-1
    } else {
      #Bin breaks
      bin.breaks<-hist(c(x,y), breaks=bw, plot=F)$breaks
      #Bin width
      band.width<-diff(bin.breaks)[1]
      #Number of bins
      bin.n<-bw
    }

    #Counting the number of elements per bin
    histx<-hist(x, breaks=bin.breaks, plot=FALSE)[[2]]
    histy<-hist(y, breaks=bin.breaks, plot=FALSE)[[2]]
    #Relative counts
    rel.histx<-histx/sum(histx)
    rel.histy<-histy/sum(histy)

    #Calculating the Bhattacharyya Coefficient (sum of the square root of the multiple of the relative counts of both distributions)
    bhatt.coeff<-sum(sqrt(rel.histx*rel.histy))
    return(bhatt.coeff)
    #End
  }


  getBc <- function(x, yhat.gene) {
    return(getBCoef(x = yhat.gene[,1,x],
                    y = yhat.gene[,2,x]))
  }


  bc <- sapply(X = 1:usage.data$N_gene, FUN = getBc,
               yhat.gene = glm.ext$Yhat_gene)

  return (bc)
}





getGroupStats <- function(glm.ext,
                          usage.data,
                          hdi.level) {


  getGroupYhat <- function(x,
                           gene.names,
                           conditions,
                           yhat.gene,
                           hdi.level,
                           usage.data) {
    hdi.1 <- getHdi(vec = yhat.gene[,1,x], hdi.level = hdi.level)
    hdi.2 <- getHdi(vec = yhat.gene[,2,x], hdi.level = hdi.level)




    # get mean raw data
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


    return(rbind(data.frame(gene_name = gene.names[x],
                            ppc.M = mean(yhat.gene[,1,x]),
                            ppc.L = hdi.1[1],
                            ppc.H = hdi.1[2],
                            raw.M = mean(real.pct.1),
                            condition = conditions[1],
                            stringsAsFactors = FALSE),
                 data.frame(gene_name = gene.names[x],
                            ppc.M = mean(yhat.gene[,2,x]),
                            ppc.L = hdi.2[1],
                            ppc.H = hdi.2[2],
                            raw.M = mean(real.pct.2),
                            condition = conditions[2],
                            stringsAsFactors = FALSE)))
  }



  getPointData <- function(usage.data) {
    pct.usage.data <- c()
    for(i in 1:ncol(usage.data$Y)) {
      row <- data.frame(gene_usage_pct = usage.data$Y[, i]/usage.data$N[i]*100,
                        condition = usage.data$Xorg[i],
                        sample_id = usage.data$sample_names[i],
                        gene_name = usage.data$gene_names,
                        stringsAsFactors = FALSE)

      pct.usage.data <- rbind(pct.usage.data, row)
    }
    return (pct.usage.data)
  }




  conditions = c(unique(usage.data$Xorg[usage.data$X == 1]),
                 unique(usage.data$Xorg[usage.data$X == -1]))

  group.ppc <- lapply(X = 1:usage.data$N_gene,
                      FUN = getGroupYhat,
                      gene.names = usage.data$gene_names,
                      condition = conditions,
                      yhat.gene = glm.ext$Yhat_gene,
                      hdi.level = hdi.level,
                      usage.data = usage.data)
  group.ppc <- do.call(rbind, group.ppc)


  individual.pct.data <- getPointData(usage.data = usage.data)



  return (list(group.ppc = group.ppc,
               individual.pct.data = individual.pct.data))
}




# two sided t.test
getTTestStats <- function(usage.data) {
  getTTest <- function(x, Ys, Xs, Ns) {
    return(try(t.test(Ys[x, ]/Ns*100~Xs)))
  }
  getTTestSummary <- function(x) {
    if(class(x) == "try-error") {
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
  tout.summary$t.test.fdr.pvalue <- p.adjust(p = tout.summary$t.test.pvalue,
                                             method = "fdr")
  tout.summary$t.test.bonf.pvalue <- p.adjust(p = tout.summary$t.test.pvalue,
                                              method = "bonferroni")

  return (tout.summary)
}





getManUStats <- function(usage.data) {

  getMTest <- function(x, Ys, Xs, Ns) {
    return(try(wilcox.test(Ys[x, ]/Ns*100~Xs)))
  }

  getMSummary <- function(x) {
    if(class(x) == "try-error") {
      return(data.frame(u.test.pvalue = NA,
                        u.test.wvalue = NA,
                        stringsAsFactors = FALSE))
    }
    return(data.frame(u.test.pvalue = x$p.value,
                      u.test.wvalue = x$statistic,
                      stringsAsFactors = FALSE))
  }

  mout <- suppressWarnings(expr = lapply(X = 1:usage.data$N_gene,
                                         FUN = getMTest,
                                         Ys = usage.data$Y,
                                         Xs = usage.data$X,
                                         Ns = usage.data$N))

  mout.summary <- do.call(rbind, lapply(mout, getMSummary))
  mout.summary$gene_name <- usage.data$gene_names

  # multiple correction
  mout.summary$u.test.fdr.pvalue <- p.adjust(p = mout.summary$u.test.pvalue,
                                             method = "fdr")
  mout.summary$u.test.bonf.pvalue <- p.adjust(p = mout.summary$u.test.pvalue,
                                              method = "bonferroni")

  return (mout.summary)
}




