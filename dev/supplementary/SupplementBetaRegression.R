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
stats$pmax.delta <- abs(stats$pmax-stats$beta_reg_pmax)

# g <- ggplot()+
#   geom_abline(slope = 1, intercept = 0, linetype = "dashed", col = "darkgray")+
#   geom_point(data = stats, aes(x = pmax, y = beta_reg_pmax), col = "red")+
#   geom_text_repel(data = stats[stats$rank.zibb <= 5 | stats$rank.zibr <= 5 |
#                                  stats$gene_name %in% c("IGHV1-58", "IGHV3-72"), ],
#                   min.segment.length = 0.1, size = 2, segment.size = 0.25,
#                   aes(x = pmax, y = beta_reg_pmax, label = gene_name_fig))+
#   theme_bw()+
#   xlab(label = expression(pi["ZIBB"]))+
#   ylab(label = expression(pi["Beta"]))
#
# g
#
# ggsave(filename = "dev/supplementary/BetaRegression.pdf", plot = g,
#        device = "pdf", width = 4, height = 4, dpi = 600)
# ggsave(filename = "dev/supplementary/BetaRegression.eps", plot = g,
#        device = "eps", width = 4, height = 4, dpi = 600)




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

s.prop$gene_name <- Mzibb$usage.data$gene_names[s.prop$g]



ggplot(data = s.prop)+
  facet_wrap(facets = ~s)+
  geom_abline(slope = 1, intercept = 0, col = "darkgray", linetype = "dashed")+
  geom_errorbar(aes(x  = zibb.M, ymin = zibb.L, ymax = zibb.H), col = "gray")+
  geom_errorbarh(aes(y = beta.M, xmin = beta.L, xmax = beta.H), col = "gray")+
  geom_point(aes(y = mean, x  = zibb.M, fill = N/10^3), shape = 21, stroke = 0.01)+
  coord_cartesian(xlim = c(-2, 2), ylim = c(-2, 2))+
  theme_bw(base_size = 9)+
  xlab(label = "ZIBB")+
  ylab(label = "Beta")



ggplot(data = s.prop[s.prop$gene_name == "IGHV3-64D", ])+
  facet_wrap(facets = ~s)+
  geom_abline(slope = 1, intercept = 0, col = "darkgray", linetype = "dashed")+
  geom_errorbar(aes(x  = zibb.M, ymin = zibb.L, ymax = zibb.H), col = "gray")+
  geom_errorbarh(aes(y = beta.M, xmin = beta.L, xmax = beta.H), col = "gray")+
  geom_point(aes(y = mean, x  = zibb.M, fill = N/10^3), shape = 21, stroke = 0.01)+
  coord_cartesian(xlim = c(-2, 2), ylim = c(-2, 2))+
  theme_bw(base_size = 9)+
  xlab(label = "ZIBB")+
  ylab(label = "Beta")




ggplot(data = s.prop)+
  facet_wrap(facets = ~s)+
  geom_density(aes(zibb.H-zibb.L), linetype = "solid")+
  geom_density(aes(beta.H-beta.L), linetype = "solid", col = "darkgray")+
  theme_bw(base_size = 9)+
  xlim(c(0, 2))+
  ggtitle(label = "black = ZIBB, gray = Beta")+
  xlab(label = "Length of 95% HDI")




e.zibb <- rstan::extract(object = Mzibb$glm)
e.beta <- rstan::extract(object = glm.prop)

hist(e.zibb$beta_gene[, 41])
hist(e.beta$beta_gene[, 41])


x <- IGHV_HCV[IGHV_HCV$gene_name == "IGHV3-64D", ]


data.list$Y[41, ]
plot(data.list$Y[41, ])


ppc <- Mzibb$ppc.data$ppc.repertoire
ppc <- ppc[ppc$gene_name == "IGHV3-64D", ]
ggplot(data = ppc)+
  geom_point(aes(x = sample_id, y = observed.count, col = condition))+
  theme_bw(base_size = 9)+
  geom_errorbar(aes(x = sample_id, y = observed.count, ymin = ppc.count.L, ymax = ppc.count.H))

ggplot(data = ppc)+
  geom_point(aes(x = sample_id, y = observed.pct, col = condition))+
  theme_bw(base_size = 9)+
  geom_errorbar(aes(x = sample_id, y = observed.pct, ymin = ppc.pct.L, ymax = ppc.pct.H))

out <- c()
for(i in 1:29) {
  hdi <- getHdi(vec = e.zibb$beta[,i,41], hdi.level = 0.95)
  row <- data.frame(M = mean(e.zibb$beta[,i,41]),
                    L = hdi[1], H = hdi[2])
  out <- rbind(out, row)
}
out$sample_id <- Mzibb$usage.data$sample_names
out$condition <- Mzibb$usage.data$Xorg
rm(i, hdi, row)


ggplot(data = out)+
  geom_errorbar(aes(x = sample_id, y = M, ymin = L, ymax = H, col = condition))+
  geom_point(aes(x = sample_id, y = M, col = condition))+
  theme_bw()

hist(e.zibb$beta_gene[, 41])
