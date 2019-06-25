require(rstan)
rstan_options(auto_write = TRUE)

source(file = "R/Util.R")
source(file = "R/Usage.R")




# D1
input <- get(load(file = "inst/input.RData"))
input.V <- input[which(regexpr(pattern = "IGHV", text = input$gene_name) != -1), ]
input.J <- input[which(regexpr(pattern = "IGHJ", text = input$gene_name) != -1), ]
rm(input)


usage.J <- compareUsage(usage.data = input.J,
                        mcmc.warmup = 500,
                        mcmc.steps = 1500,
                        mcmc.chains = 2,
                        mcmc.cores = 2,
                        hdi.level = 0.95,
                        adapt.delta = 0.95,
                        max.treedepth = 13)

usage.V <- compareUsage(usage.data = input.V,
                        mcmc.warmup = 500,
                        mcmc.steps = 1500,
                        mcmc.chains = 2,
                        mcmc.cores = 2,
                        hdi.level = 0.95,
                        adapt.delta = 0.95,
                        max.treedepth = 13)






# J-viz
ggplot()+
  geom_errorbar(data = usage.J$group.ppc$group.ppc,
                aes(x = gene_name, ymin = ppc.L,
                    ymax = ppc.H, col = condition),
                position = position_dodge(width = .8),
                width = 0.75)+
  geom_point(data = usage.J$group.ppc$individual.pct.data,
             aes(x = gene_name, y = gene_usage_pct,
                 col = condition), shape = 21, size = 0.5,
             position = position_dodge(width = .8), fill = "black")+
  theme_bw(base_size = 9)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4),
        legend.position = "top")







# D2
input <- get(load(file = "inst/IGHV_HCV_CSM.RData"))
input <- input[, c("sample_id", "condition", "gene_name", "gene_usage_count")]
input.V <- input


usage.V <- compareUsage(usage.data = input.V,
                        mcmc.warmup = 500,
                        mcmc.steps = 1500,
                        mcmc.chains = 4,
                        mcmc.cores = 4,
                        hdi.level = 0.95,
                        adapt.delta = 0.95,
                        max.treedepth = 13)

z <- usage.V$glm
z <- summary(z, pars = c("alpha_gene", "beta_gene"))$summary
plot(z[1:(nrow(z)/2), 1], z[(nrow(z)/2 + 1):nrow(z), 1])




z <- usage.V$glm.summary


ggplot()+
  geom_errorbar(data = usage.V$group.ppc$group.ppc,
                aes(x = gene_name, ymin = ppc.L,
                    ymax = ppc.H, col = condition),
                position = position_dodge(width = .8),
                width = 0.75)+
  geom_point(data = usage.V$group.ppc$individual.pct.data,
             aes(x = gene_name, y = gene_usage_pct,
                 col = condition), shape = 21, size = 0.5,
             position = position_dodge(width = .8), fill = "black")+
  geom_text(data = usage.V$glm.summary, size = 2,
            aes(x = gene, y = 20, label = round(x = pmax, digits = 2)))+
  theme_bw(base_size = 9)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4),
        legend.position = "top")


ppc <- usage.V$ppc
ggplot(data = ppc)+
  facet_wrap(facets = ~sample_id, ncol = 4)+
  geom_point(aes(x = raw, y = ppc.raw.mean, col = condition))+
  geom_abline(intercept = 0, slope = 1)

group.ppc <- usage.V$group.ppc$individual.pct.data
ggplot(data = ppc[ppc$sample_id == "HD1", ])+
  facet_wrap(facets = ~sample_id, ncol = 4)+
  geom_point(aes(x = raw, y = ppc.raw.mean, col = condition))+
  geom_errorbar(aes(x = raw, y = ppc.raw.mean, ymin = ppc.raw.L, ymax = ppc.raw.H, col = condition))+
  geom_abline(intercept = 0, slope = 1)


usage.V$t.test.stats
hist(usage.V$man.u.stats$fdr.p.value, breaks = 100)
sort(usage.V$man.u.stats$fdr.p.value)

#
# e <- rstan::extract(usage.V$glm, pars = c("alpha_gene", "beta_gene"))
# plot(density(1/(1 + exp(-(e$alpha_gene[,26]+e$beta_gene[,26])))))
# points(density(1/(1 + exp(-(e$alpha_gene[,26]+e$beta_gene[,26]*(-1))))),
#        col = "red", type = "l")
#
# dim(e$beta_gene)
# hist(e$beta_gene[, 26], breaks = 100)
# usage.V$glm.summary[26,]
#
# usage.V$t.test.stats
#
#
#
# z <- usage.V$glm
# z <- summary(z, pars = c("alpha_gene", "beta_gene"))$summary
# plot(z[1:65, 1], z[66:nrow(z), 1])
