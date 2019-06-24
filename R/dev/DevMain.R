require(rstan)
rstan_options(auto_write = TRUE)

source(file = "R/Util.R")
source(file = "R/Usage.R")


input <- get(load(file = "inst/input.RData"))
input.V <- input[which(regexpr(pattern = "IGHV", text = input$gene_name) != -1), ]
input.J <- input[which(regexpr(pattern = "IGHJ", text = input$gene_name) != -1), ]
rm(input)


comparison.usage <- compareUsage(usage.data = input.J,
                                 mcmc.warmup = 500,
                                 mcmc.steps = 1500,
                                 mcmc.chains = 2,
                                 mcmc.cores = 2,
                                 hdi.level = 0.95,
                                 adapt.delta = 0.95,
                                 max.treedepth = 13)

comparison.usage$glm.summary









ggplot()+
  geom_errorbar(data = comparison.usage$group.ppc$group.ppc,
                aes(x = gene_name, ymin = ppc.L,
                    ymax = ppc.H, col = condition),
                position = position_dodge(width = .8),
                width = 0.75)+
  geom_point(data = comparison.usage$group.ppc$individual.pct.data,
             aes(x = gene_name, y = gene_usage_pct,
                 col = condition), shape = 21, size = 0.5,
             position = position_dodge(width = .8), fill = "black")+
  theme_bw(base_size = 9)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4),
        legend.position = "top")


comparison.usage$usage.data$
t.test()



comparison.usage$glm.summary

