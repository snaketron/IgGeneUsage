require(ggplot2)
require(gridExtra)
require(ggrepel)
require(rstan)
rstan_options(auto_write = TRUE)

source(file = "R/Util.R")
source(file = "R/Usage.R")




IGHV_HCV <- read.csv(file = "inst/IGHV_HCV_CSM.csv", sep = " ", as.is = T)
IGHV_HCV <- IGHV_HCV[, c("sample_id", "condition",
                         "gene_name", "gene_usage_count")]


IGHV_HCV <- IGHV_HCV[which(regexpr(pattern = "IGHV3",
                                   text = IGHV_HCV$gene_name) != -1), ]


stan.files <- list.files(path = "src/stan_files/",
                         pattern = "\\.stan",
                         full.names = T)
stan.files.short <- list.files(path = "src/stan_files/",
                         pattern = "\\.stan",
                         full.names = F)

for(m in 1:length(stan.files)) {
  M <- diffUsage(usage.data = IGHV_HCV,
                             mcmc.warmup = 1000,
                             mcmc.steps = 2500,
                             mcmc.chains = 4,
                             mcmc.cores = 4,
                             hdi.level = 0.95,
                             adapt.delta = 0.95,
                             max.treedepth = 13,
                             dev.model = stan.files[m])

  out <- paste("R/dev/ighv_hcv_", gsub(pattern = "\\.stan",
                                       replacement = '\\.RData',
                                       x = stan.files.short[m]), sep = '')
  save(M, file = out)
}
