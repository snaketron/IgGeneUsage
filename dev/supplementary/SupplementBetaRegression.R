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
data.list$Yp <- Yp
rm(Yp, i)



# model.prop <- rstan::stan_model(file = "src/stan_files/case_by_case/proportion_h_multigene.stan")
# glm.prop <- rstan::sampling(object = model.prop,
#                             data = data.list,
#                             chains = 4,
#                             cores = 4,
#                             iter = 5000,
#                             warmup = 2500,
#                             control = list(adapt_delta = 0.95,
#                                            max_treedepth = 10))
# save(glm.prop, file = "dev/supplementary/proportion_h_multigene.RData")



glm.prop <- get(load(file = "dev/supplementary/proportion_h_multigene.RData"))
s <- data.frame(summary(glm.prop, par = "beta_gene")$summary)
s$gene_name <- data.list$gene_names
e <- rstan::extract(object = glm.prop)
s$pmax <- getPmax(glm.ext = e)
colnames(s) <- paste("beta_reg", colnames(s), sep = '_')




# zibb
Mzibb <- get(load(file = "dev/ighv_hcv_zibb_model.RData"))
z <- Mzibb$glm.summary
d <- merge(x = z, y = s, by.x = "gene_name", by.y = "beta_reg_gene_name")



# format data
stats <- Mzibb$test.stats
stats <- merge(x = d, y = stats, by = "gene_name")


# ranks
stats <- stats[order(stats$pmax, decreasing = T), ]
stats$rank.zibb <- 1:nrow(stats)
stats <- stats[order(stats$beta_reg_pmax, decreasing = T), ]
stats$rank.zibr <- 1:nrow(stats)
stats <- stats[order(stats$t.test.fdr.pvalue, decreasing = F), ]
stats$rank.t <- 1:nrow(stats)
stats <- stats[order(stats$u.test.fdr.pvalue, decreasing = F), ]
stats$rank.u <- 1:nrow(stats)
stats$gene_name_fig <- gsub(pattern = "IGHV", replacement = '', x = stats$gene_name)


g <- ggplot()+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", col = "darkgray")+
  geom_point(data = stats, aes(x = pmax, y = beta_reg_pmax), col = "red")+
  geom_text_repel(data = stats[stats$rank.zibb <= 5 | stats$rank.zibr <= 5 |
                                 stats$gene_name %in% c("IGHV1-58", "IGHV3-72"), ],
                  min.segment.length = 0.1, size = 2, segment.size = 0.25,
                  aes(x = pmax, y = beta_reg_pmax, label = gene_name_fig))+
  theme_bw()+
  xlab(label = expression(pi["ZIBB"]))+
  ylab(label = expression(pi["Beta"]))

g



ggsave(filename = "dev/supplementary/BetaRegression.pdf", plot = g,
       device = "pdf", width = 4, height = 4, dpi = 600)
ggsave(filename = "dev/supplementary/BetaRegression.eps", plot = g,
       device = "eps", width = 4, height = 4, dpi = 600)







# effect on individual betas
s.zibb <- summary(Mzibb$glm, pars = "beta")$summary
s.prop <- summary(glm.prop, pars = "beta_sample")$summary

s.prop <- data.frame(s.prop)
s.prop$zibb.M <- s.zibb[, 1]
s.prop$zibb.L <- s.zibb[, 4]
s.prop$zibb.H <- s.zibb[, 8]

s.prop$beta.M <- s.prop$mean
s.prop$beta.L <- s.prop$X2.5.
s.prop$beta.H <- s.prop$X97.5.

s.prop$g <- rep(x = 1:69, times = 29)
s.prop$s <- rep(1:29, each  = 69)
s.prop$N <- NA
for(i in 1:29) {
  s.prop$N[s.prop$s == i] <- Mzibb$usage.data$N[i]
}

ggplot(data = s.prop)+
  facet_wrap(facets = ~s)+
  geom_abline(slope = 1, intercept = 0, col = "darkgray", linetype = "dashed")+
  geom_errorbar(aes(y = beta.M, x  = zibb.M, ymin = zibb.L, ymax = zibb.H), col = "gray")+
  geom_errorbarh(aes(y = beta.M, x  = zibb.M, xmin = beta.L, xmax = beta.H), col = "gray")+
  geom_point(aes(y = mean, x  = zibb.M, fill = N/10^3), shape = 21, stroke = 0.01)+
  coord_cartesian(xlim = c(-2, 2), ylim = c(-2, 2))+
  theme_bw()

ggplot(data = s.prop)+
  facet_wrap(facets = ~s)+
  geom_density(aes(zibb.H-zibb.L), linetype = "solid")+
  geom_density(aes(beta.H-beta.L), linetype = "solid", col = "darkgray")+
  theme_bw()+
  xlim(c(0, 2))+
  ggtitle(label = "black = ZIBB, gray = Beta")
