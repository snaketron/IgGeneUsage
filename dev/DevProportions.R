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
# Yp[which(Yp == 0, arr.ind = T)] <- 10^(-6) # numerical stability log(0)
data.list$Yp <- Yp
rm(Yp, i)


# model.proportions <- rstan::stan_model(file = "src/stan_files/case_by_case/proportion_multigene.stan")
model.h.proportions <- rstan::stan_model(file = "src/stan_files/case_by_case/proportion_h_multigene.stan")


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
e <- rstan::extract(object = glm.h.prop)
s$pmax <- getPmax(glm.ext = e)
hist(s$mean, breaks = 100)
hist(s$pmax, breaks = 100)
colnames(s) <- paste("beta.reg", colnames(s), sep = '_')


# zibb
Mzibb <- get(load(file = "dev/ighv_hcv_zibb_model.RData"))
z <- Mzibb$glm.summary

d <- merge(x = z, y = s, by.x = "gene_name", by.y = "beta.reg_gene_name")
d$delta <- d$pmax-d$beta.reg_pmax
d$delta_effect <- abs(d$effect_mean)-abs(d$beta.reg_mean)

ggplot()+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", col = "darkgray")+
  geom_point(data = d, aes(x = pmax, y = beta.reg_pmax))+
  geom_text_repel(data = d[d$pmax >= 0.8 | d$beta.reg_pmax >= 0.8, ],
                  aes(x = pmax, y = beta.reg_pmax, label = gene_name))+
  theme_bw()

ggplot()+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", col = "darkgray")+
  geom_point(data = d, aes(x = abs(effect_mean), y = abs(beta.reg_mean)))+
  geom_errorbar(data = d, aes(x = abs(effect_mean), y = abs(beta.reg_mean),
                              ymin = beta.reg_X2.5., ymax = beta.reg_X97.5.))+
  geom_text_repel(data = d[abs(d$delta_effect) >= 0.02, ],
                  aes(x = abs(effect_mean), y = abs(beta.reg_mean), label = gene_name))+
  theme_bw()
