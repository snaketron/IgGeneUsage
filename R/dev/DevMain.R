require(rstan)
rstan_options(auto_write = TRUE)

source(file = "R/Util.R")
source(file = "R/Usage.R")




# D1
input <- get(load(file = "inst/input.RData"))
input.V <- input[which(regexpr(pattern = "IGHV", text = input$gene_name) != -1), ]
input.J <- input[which(regexpr(pattern = "IGHJ", text = input$gene_name) != -1), ]
rm(input)


# D2
input <- get(load(file = "inst/IGHV_HCV.RData"))
input <- input[, c("sample_id", "condition", "gene_name", "gene_usage_count")]
input.V <- input


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
                        max.treedepth = 10)









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
  theme_bw(base_size = 9)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4),
        legend.position = "top")

usage.V$t.test.stats
usage.V$man.u.stats



z <- usage.V$glm
z <- summary(z, pars = c("alpha_gene", "beta_gene"))$summary
plot(z[1:65, 1], z[66:nrow(z), 1])
