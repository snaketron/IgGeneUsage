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
Yp[which(Yp == 0, arr.ind = T)] <- 10^(-16) # numerical stability log(0)
data.list$Yp <- Yp
rm(Yp, i)


model.proportions <- rstan::stan_model(file = "src/stan_files/case_by_case/proportion_multigene.stan")
model.h.proportions <- rstan::stan_model(file = "src/stan_files/case_by_case/proportion_h_multigene.stan")



glm.prop <- rstan::sampling(object = model.proportions,
                            data = data.list,
                            chains = 4,
                            cores = 4,
                            iter = 3000,
                            warmup = 1000,
                            control = list(adapt_delta = 0.95,
                                           max_treedepth = 10))


glm.h.prop <- rstan::sampling(object = model.h.proportions,
                            data = data.list,
                            chains = 4,
                            cores = 4,
                            iter = 3000,
                            warmup = 1000,
                            control = list(adapt_delta = 0.95,
                                           max_treedepth = 10))


s <- data.frame(summary(glm.h.prop, par = "beta_gene")$summary)
s$gene_name <- data.list$gene_names


Mzibb <- get(load(file = "dev/ighv_hcv_zibb_model.RData"))
z <- Mzibb$glm.summary
z$gene_name == s$gene_name

z$p.mean <- s$mean
z$p.median <- s$X50.
e <- rstan::extract(object = glm.h.prop)
z$p.pmax <- getPmax(glm.ext = e)
z$delta <- z$pmax - z$p.pmax

plot(abs(z$effect_mean), abs(z$p))


plot(abs(z$pmax), abs(z$p.pmax))
