model <- get(load(file = "R/dev/model.RData"))
zib.s <- get(load(file = "R/dev/zib_singlez.RData"))
zib.m <- get(load(file = "R/dev/zib_multiz.RData"))
beta.bin <- get(load(file = "R/dev/beta_binomial.RData"))
rm(usage.hcv.csm)


loo::compare(loo::loo(loo::extract_log_lik(model$glm)),
             loo::loo(loo::extract_log_lik(zib.s$glm)),
             loo::loo(loo::extract_log_lik(zib.m$glm)),
             loo::loo(loo::extract_log_lik(beta.bin$glm)))


plot(abs(b$glm.summary[, 1]), abs(zib.m$glm.summary[, 1]))
abline(0, 1)

which(abs(b$glm.summary[, 1]) > 0.6 & abs(zib.m$glm.summary[, 1]) < 0.1)
gene_name <- b$glm.summary$gene_name[65]


# gene_name <- "IGHV3-30-3"
group.ppc <- b$group.ppc$group.ppc
individual.data <- b$group.ppc$individual.pct.data

group.ppc <- zib.m$group.ppc$group.ppc
individual.data <- zib.m$group.ppc$individual.pct.data


group.ppc <- group.ppc[group.ppc$gene_name == gene_name, ]
individual.data <- individual.data[individual.data$gene_name == gene_name, ]

ggplot()+
  geom_errorbar(data = group.ppc,
                aes(x = gene_name, ymin = ppc.L, ymax = ppc.H, col = condition),
                position = position_dodge(width = .8), width = 0.75)+
  geom_point(data = individual.data,
             aes(x = gene_name, y = gene_usage_pct, col = condition),
             shape = 21, size = 2, fill = "black",
             position = position_jitterdodge(jitter.width = 0.15,
                                             jitter.height = 0,
                                             dodge.width = 0.8))+
  theme_bw(base_size = 9)+
  theme(legend.position = "top")
