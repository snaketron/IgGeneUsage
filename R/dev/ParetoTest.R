require(ggplot2)
require(gridExtra)
require(ggrepel)
require(rstan)
rstan_options(auto_write = TRUE)

source(file = "R/Util.R")
source(file = "R/Usage.R")



J <- get(load(file = "inst/temp_usage.RData"))
rm(usage, nb.usage)
J <- J[which(regexpr(pattern = "IGHJ", text = J$Seg) != -1), ]
J$sample_id <- J$sample
J$condition <- J$key
J$gene_name <- J$Seg
J$gene_usage_count <- J$cell.count
IGHJ_HCV <- J
rm(J)

IGHJ_HCV <- IGHJ_HCV[, c("sample_id", "condition",
                         "gene_name", "gene_usage_count")]

M <- diffUsage(usage.data = IGHJ_HCV,
               mcmc.warmup = 1000,
               mcmc.steps = 2500,
               mcmc.chains = 4,
               mcmc.cores = 4,
               hdi.level = 0.95,
               adapt.delta = 0.95,
               max.treedepth = 13,
               dev.model = "src/stan_files/beta_binomial_model.stan")


M <- get(load(file = "R/dev/ighv_hcv_beta_binomial_model.RData"))
loo.M <- loo::loo(loo::extract_log_lik(M$glm))
loo.M
plot(loo.M$diagnostics$pareto_k)

L <- matrix(data = F, nrow = nrow(M$usage.data$Y), ncol = ncol(M$usage.data$Y))
L[which(loo.M$diagnostics$pareto_k >= 1)] <- T

q <- c()
for(i in 1:nrow(L)) { # i = gene
  for(j in 1:ncol(L)) { # j = sample
    if(L[i,j]) {
      row <- data.frame(Y = M$usage.data$Y[i,j],
                        N = M$usage.data$N[j],
                        pct = M$usage.data$Y[i,j]/M$usage.data$N[j]*100,
                        gene_name = M$usage.data$gene_names[i],
                        sample_id = M$usage.data$sample_names[j])
      q <- rbind(q, row)
    }
  }
}
rm(L, i, j, row)

ggplot()+
  geom_errorbar(data = M$group.ppc$group.ppc,
                aes(x = gene_name, ymin = ppc.L,
                    ymax = ppc.H, col = condition),
                position = position_dodge(width = .8), width = 0.75)+
  geom_point(data = M$group.ppc$individual.pct.data,
             aes(x = gene_name, y = gene_usage_pct, col = condition),
             shape = 21, size = 1.5, fill = "black",
             position = position_jitterdodge(jitter.width = 0.25,
                                             jitter.height = 0,
                                             dodge.width = 0.8))+
  geom_point(data = q, aes(x = gene_name, y = pct), col = "blue", size = 1.5)+
  theme_bw(base_size = 9)+
  theme(legend.position = "top")


ilogit <- function(x) {return(1/(1+exp(-(x))))}
b <- function(x,n) {rbinom(n = 1, size = n, prob = x)}
e <- rstan::extract(object = M$glm)

g <- 26
x <- (-1)
n <- median(M$usage.data$N[M$usage.data$X == x])
p <- sapply(ilogit(e$alpha_gene[, g]+e$beta_gene[, g]*x), b, n = n)/n * 100
hist(p)
getHdi(vec = p, hdi.level = 0.99)
