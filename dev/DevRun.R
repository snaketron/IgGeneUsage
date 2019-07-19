require(ggplot2)
require(gridExtra)
require(ggrepel)
require(rstan)
rstan_options(auto_write = TRUE)

source(file = "R/Util.R")
source(file = "R/Usage.R")


# HCV
data("IGHV_HCV")


# # HIV
# data("IGHV_HIV")
#
#
# # Alakazam
# data("Ig")
#
#
# # Epitopes
# data("CDR3_Epitopes")





d <- IGHV_HCV
# d <- CDR3_Epitopes




stan.files <- list.files(path = "src/stan_files/",
                         pattern = "\\.stan",
                         full.names = T)
stan.files.short <- list.files(path = "src/stan_files/",
                               pattern = "\\.stan",
                               full.names = F)


for(m in 1:length(stan.files)) {
  M <- diffUsage(usage.data = d,
                 mcmc.warmup = 500,
                 mcmc.steps = 1500,
                 mcmc.chains = 4,
                 mcmc.cores = 4,
                 hdi.level = 0.95,
                 adapt.delta = 0.95,
                 max.treedepth = 13,
                 dev.model = stan.files[m])

  # out <- paste("dev/epitopes_", gsub(pattern = "\\.stan",
  #                                      replacement = '\\.RData',
  #                                      x = stan.files.short[m]), sep = '')
  # out <- paste("dev/ighv_hiv_", gsub(pattern = "\\.stan",
  #                                      replacement = '\\.RData',
  #                                      x = stan.files.short[m]), sep = '')
  out <- paste("dev/ighv_hcv_", gsub(pattern = "\\.stan",
                                     replacement = '\\.RData',
                                     x = stan.files.short[m]), sep = '')
  # out <- paste("dev/alakazam_", gsub(pattern = "\\.stan",
  #                                      replacement = '\\.RData',
  #                                      x = stan.files.short[m]), sep = '')
  save(M, file = out)
}
