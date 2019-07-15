require(ggplot2)
require(gridExtra)
require(ggrepel)
require(rstan)
rstan_options(auto_write = TRUE)

source(file = "R/Util.R")
source(file = "R/Usage.R")


# HCV
data("IGHV_HCV")


data.list <- getUsageData(usage = IGHV_HCV)
Yp <- data.list$Y
for(i in 1:nrow(Yp)) {
  Yp[i, ] <- Yp[i, ]/as.numeric(data.list$N)
}
Yp[which(Yp == 0, arr.ind = T), ] <- 10^(-16) # numerical stability log(0)
data.list$Yp <- Yp
rm(Yp, i)


model.proportions <- rstan::stan_model(file = "src/stan_files/case_by_case/proportion_multigene.stan")



glm.prop <- rstan::sampling(object = model.proportions,
                            data = data.list,
                            chains = 4,
                            cores = 4,
                            iter = 1000,
                            warmup = 500,
                            control = list(adapt_delta = 0.90,
                                           max_treedepth = 10))


summary(glm.prop, par = "beta_gene")$summary




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
