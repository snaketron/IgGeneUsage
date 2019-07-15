require(ggplot2)
require(gridExtra)
require(ggrepel)
require(rstan)
rstan_options(auto_write = TRUE)

source(file = "R/Util.R")
source(file = "R/Usage.R")


# HCV
data("IGHV_HCV")
# d <- IGHV_HCV[IGHV_HCV$gene_name %in% c("IGHV1-58"), ]
# d <- IGHV_HCV[IGHV_HCV$gene_name %in% c("IGHV1-58", "IGHV3-72"), ]


model.binomial <- rstan::stan_model(file = "src/stan_files/case_by_case//binomial.stan")
model.proportions <- rstan::stan_model(file = "src/stan_files/case_by_case/proportion.stan")
model.beta.binomial <- rstan::stan_model(file = "src/stan_files/case_by_case/beta_binomial.stan")


data.list <- getUsageData(usage = IGHV_HCV)
i <- which(data.list$gene_names == "IGHV1-58")
# i <- which(data.list$gene_names == "IGHV3-72")
data.list$Y <- data.list$Y[i, ]
data.list$gene_names <- data.list$gene_names[i]

glm <- rstan::sampling(object = bm,
                       data = data.list,
                       chains = 4,
                       cores = 4,
                       iter = 5000,
                       warmup = 1500,
                       refresh = 500,
                       control = list(adapt_delta = 0.999))

summary(glm, par = "b_sample")$summary
summary(glm, par = "b")$summary

e <- rstan::extract(object = glm)

par(mfrow = c(2, 2))
hist(e$Yhat_group[, 1]*100)
hist(e$Yhat_group[, 2]*100)


plot(data.list$Y/data.list$N*100, col = as.factor(data.list$X))
hist(e$Yhat_individual[, 15])

hist(e$b_sample, breaks = 100)
hist(e$Yhat_individual[, 5], breaks = 100)

sum(e$b_sample < 0)/length(e$b_sample)


plot(density(100*1/(1 + exp(-(e$a_sample + e$b_sample*1)))), xlim = c(0, 0.5))
points(density(100*1/(1 + exp(-(e$a_sample + e$b_sample*(-1))))), xlim = c(0, 0.5), type = "l")
plot(d$gene_usage_count, col = as.factor(d$condition))
data.list$N


getHdi(vec = e$Yhat[, u], hdi.level = 0.95)
