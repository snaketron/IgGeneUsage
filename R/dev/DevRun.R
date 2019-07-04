require(ggplot2)
require(gridExtra)
require(ggrepel)
require(rstan)
rstan_options(auto_write = TRUE)

source(file = "R/Util.R")
source(file = "R/Usage.R")


# HCV
# d <- read.csv(file = "inst/IGHV_HCV_CSM.csv", sep = " ", as.is = T)
# d <- d[, c("sample_id", "condition",
#            "gene_name", "gene_usage_count")]

# HIV
# d <- read.csv(file = "inst/IGHV_HIV.csv", sep = " ", as.is = T)
# d <- d[, c("sample_id", "condition",
#            "gene_name", "gene_usage_count")]


# TRV_Cancer -> some processing needed before use, too unreliable
# d <- read.csv(file = "inst/TCR_Cancer.csv", sep = ";", as.is = T)
# d <- d[d$GroupInformation %in% c("Tumors", "None Tumor"), ]
# d <- d[which(regexpr(pattern = "TRBV", text = d$ReferenceName) != -1),]
# d$sample_id <- paste(d$IndividualID, d$GroupInformation, sep = '_')
# d$condition <- d$GroupInformation
# d$gene_name <- d$ReferenceName
# d$gene_usage_count <- d$UsageNumber
# d <- d[, c("sample_id", "condition",
#            "gene_name", "gene_usage_count")]


# MS
# d <- read.csv(file = "inst/TRBV_MS.csv", sep = " ", as.is = T)
# d$gene_usage_count <- as.integer(d$gene_usage_count)



# Ig from alakazam
d <- read.csv(file = "inst/Ig.csv", sep = " ", as.is = T)
d$gene_usage_count <- as.integer(d$gene_usage_count)




stan.files <- list.files(path = "src/stan_files/",
                         pattern = "\\.stan",
                         full.names = T)
stan.files.short <- list.files(path = "src/stan_files/",
                               pattern = "\\.stan",
                               full.names = F)


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

  out <- paste("R/dev/alakazam_ighv_families_", gsub(pattern = "\\.stan",
                                         replacement = '\\.RData',
                                         x = stan.files.short[m]), sep = '')
  save(M, file = out)
}
