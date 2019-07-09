require(ggplot2)
require(gridExtra)
require(ggrepel)
require(rstan)
rstan_options(auto_write = TRUE)

source(file = "R/Util.R")
source(file = "R/Usage.R")


# HCV
d <- read.csv(file = "inst/IGHV_HCV_CSM.csv", sep = " ", as.is = T)
d <- d[, c("sample_id", "condition",
           "gene_name", "gene_usage_count")]

# HIV
# d <- read.csv(file = "inst/IGHV_HIV.csv", sep = " ", as.is = T)
# d <- d[, c("sample_id", "condition",
#            "gene_name", "gene_usage_count")]



# Ig from alakazam
# d <- read.csv(file = "inst/Ig.csv", sep = " ", as.is = T)
# d$gene_usage_count <- as.integer(d$gene_usage_count)




stan.files <- list.files(path = "src/stan_files/",
                         pattern = "\\.stan",
                         full.names = T)
stan.files.short <- list.files(path = "src/stan_files/",
                               pattern = "\\.stan",
                               full.names = F)
m <- 6

for(m in 1:length(stan.files)) {
  M <- diffUsage(usage.data = d,
                 mcmc.warmup = 1000,
                 mcmc.steps = 2500,
                 mcmc.chains = 4,
                 mcmc.cores = 4,
                 hdi.level = 0.95,
                 adapt.delta = 0.95,
                 max.treedepth = 13,
                 dev.model = stan.files[m])

  # out <- paste("R/dev/alakazam_ighv_families_", gsub(pattern = "\\.stan",
  #                                        replacement = '\\.RData',
  #                                        x = stan.files.short[m]), sep = '')
  out <- paste("R/dev/ighv_hcv_zmix_", gsub(pattern = "\\.stan",
                                            replacement = '\\.RData',
                                            x = stan.files.short[m]), sep = '')
  save(M, file = out)
}
