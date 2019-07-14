require(ggplot2)
require(gridExtra)
require(ggrepel)
require(rstan)
rstan_options(auto_write = TRUE)

source(file = "R/Util.R")
source(file = "R/Usage.R")


# HCV
data("IGHV_HCV")



# Ig from alakazam
# data("Ig")
# d$gene_usage_count <- as.integer(d$gene_usage_count)



Mbb <- diffUsage(usage.data = IGHV_HCV,
                 mcmc.warmup = 1000,
                 mcmc.steps = 2500,
                 mcmc.chains = 4,
                 mcmc.cores = 4,
                 hdi.level = 0.95,
                 adapt.delta = 0.99,
                 max.treedepth = 10,
                 dev.model = "src/stan_files/simple_models/simple_beta_binomial_model.stan")





Mzibb <- diffUsage(usage.data = IGHV_HCV,
                   mcmc.warmup = 1000,
                   mcmc.steps = 2500,
                   mcmc.chains = 4,
                   mcmc.cores = 4,
                   hdi.level = 0.95,
                   adapt.delta = 0.99,
                   max.treedepth = 10,
                   dev.model = "src/stan_files/simple_models/simple_zibb_model.stan")



Mzibb2 <- diffUsage(usage.data = IGHV_HCV,
                   mcmc.warmup = 1000,
                   mcmc.steps = 2500,
                   mcmc.chains = 4,
                   mcmc.cores = 4,
                   hdi.level = 0.95,
                   adapt.delta = 0.99,
                   max.treedepth = 10,
                   dev.model = "src/stan_files/simple_models/simple_zibb2_model.stan")



Mmix <- diffUsage(usage.data = IGHV_HCV,
                 mcmc.warmup = 1000,
                 mcmc.steps = 2500,
                 mcmc.chains = 4,
                 mcmc.cores = 4,
                 hdi.level = 0.95,
                 adapt.delta = 0.99,
                 max.treedepth = 10,
                 dev.model = "src/stan_files/simple_models/simple_zmix_model.stan")



loo::loo(Mbb$glm)

loo::loo(Mzibb$glm)

# loo::loo(Mzibb2$glm)

loo::loo(Mmix$glm)

loo::compare(loo::loo(Mbb$glm),
             loo::loo(Mzibb$glm),
             loo::loo(Mzibb2$glm),
             loo::loo(Mmix$glm))



# 1. VIZ data
# we can compute the total number of rearrangements per sample
total.usage <- aggregate(gene_usage_count~sample_id+condition,
                         FUN = sum, data = IGHV_HCV)
total.usage$total <- total.usage$gene_usage_count
total.usage$gene_usage_count <- NULL

# merge it with the original data
viz <- merge(x = IGHV_HCV, y = total.usage,
             by = c("sample_id", "condition"),
             all.x = TRUE)
# compute %
viz$gene_usage_pct <- viz$gene_usage_count/viz$total*100






#### 2. Model parameters ####
# Mmix
Mmix.stats <- Mmix$glm.summary
Mmix.stats$effect_col <- ifelse(test = Mmix.stats$effect_L <= 0 & Mmix.stats$effect_H >= 0,
                                yes = "include", no = "exclude")
Mmix.stats <- Mmix.stats[order(abs(Mmix.stats$effect_mean), decreasing = F), ]
Mmix.stats$gene_fac <- factor(x = Mmix.stats$gene_name,
                              levels = Mmix.stats$gene_name)
Mmix.stats <- merge(x = Mmix.stats, y = Mmix$test.stats, by = "gene_name")
Mmix.stats$model <- "Mix"


# Mzibb
Mzibb.stats <- Mzibb$glm.summary
Mzibb.stats$effect_col <- ifelse(test = Mzibb.stats$effect_L <= 0 & Mzibb.stats$effect_H >= 0,
                               yes = "include", no = "exclude")
Mzibb.stats <- Mzibb.stats[order(abs(Mzibb.stats$effect_mean), decreasing = F), ]
Mzibb.stats$gene_fac <- factor(x = Mzibb.stats$gene_name,
                             levels = Mzibb.stats$gene_name)
Mzibb.stats <- merge(x = Mzibb.stats, y = Mzibb$test.stats, by = "gene_name")
Mzibb.stats$model <- "ZIBB"


# Mzibb2
Mzibb2.stats <- Mzibb2$glm.summary
Mzibb2.stats$effect_col <- ifelse(test = Mzibb2.stats$effect_L <= 0 & Mzibb2.stats$effect_H >= 0,
                                 yes = "include", no = "exclude")
Mzibb2.stats <- Mzibb2.stats[order(abs(Mzibb2.stats$effect_mean), decreasing = F), ]
Mzibb2.stats$gene_fac <- factor(x = Mzibb2.stats$gene_name,
                               levels = Mzibb2.stats$gene_name)
Mzibb2.stats <- merge(x = Mzibb2.stats, y = Mzibb2$test.stats, by = "gene_name")
Mzibb2.stats$model <- "ZIBB-2"


# Mbb
Mbb.stats <- Mbb$glm.summary
Mbb.stats$effect_col <- ifelse(test = Mbb.stats$effect_L <= 0 & Mbb.stats$effect_H >= 0,
                               yes = "include", no = "exclude")
Mbb.stats <- Mbb.stats[order(abs(Mbb.stats$effect_mean), decreasing = F), ]
Mbb.stats$gene_fac <- factor(x = Mbb.stats$gene_name,
                             levels = Mbb.stats$gene_name)
Mbb.stats <- merge(x = Mbb.stats, y = Mbb$test.stats, by = "gene_name")
Mbb.stats$model <- "Beta-binomial"


stats <- rbind(Mmix.stats, Mzibb2.stats, Mzibb.stats, Mbb.stats)
rm(Mmix.stats, Mzibb2.stats, Mzibb.stats, Mbb.stats)

# visualize disruption effects on each gene
ggplot(data = stats)+
  geom_vline(xintercept = 0, linetype = "dashed", col = "darkgray")+
  geom_errorbarh(aes(y = gene_fac, xmin = effect_L,
                     xmax = effect_H, col = model))+
  geom_point(aes(y = gene_fac, x = effect_mean, col = model))+
  theme_bw(base_size = 9)+
  theme(legend.position = "top")+
  xlab(label = "Mean effect")



# now select interesting genes
promising.genes <- Mmix$usage.data$gene_names








#### 3. Model predictions (individual level) ####
Mmix.individual <- Mmix$ppc
Mmix.individual$model <- "Mix"

Mzibb.individual <- Mzibb$ppc
Mzibb.individual$model <- "ZIBB"

Mzibb2.individual <- Mzibb2$ppc
Mzibb2.individual$model <- "ZIBB-2"

Mbb.individual <- Mbb$ppc
Mbb.individual$model <- "Beta-binomial"

individual <- rbind(Mmix.individual, Mzibb.individual, Mzibb2.individual, Mbb.individual)
rm(Mmix.individual, Mzibb.individual, Mzibb2.individual, Mbb.individual)

ggplot()+
  facet_grid(facets = model~gene_name, scales = "free_x")+
  geom_errorbar(data = individual[individual$gene_name %in% promising.genes, ],
                aes(x = sample_id, ymin = ppc.pct.L,
                    ymax = ppc.pct.H, col = condition),
                position = position_dodge(width = .8), width = 0.75)+
  geom_point(data = viz[viz$gene_name %in% promising.genes, ],
             aes(x = sample_id, y = gene_usage_pct, fill = condition,
                 shape = condition))+
  theme_bw(base_size = 9)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4))+
  theme(legend.position = "top")





#### 4. Model predictions (gene level) ####
Mmix.group <- Mmix$group.ppc$group.ppc
Mmix.group$model <- "Mix"

Mzibb.group <- Mzibb$group.ppc$group.ppc
Mzibb.group$model <- "ZIBB"

Mzibb2.group <- Mzibb2$group.ppc$group.ppc
Mzibb2.group$model <- "ZIBB-2"

Mbb.group <- Mbb$group.ppc$group.ppc
Mbb.group$model <- "Beta-binomial"

group <- rbind(Mmix.group, Mzibb.group, Mzibb2.group, Mbb.group)
rm(Mmix.group, Mzibb.group, Mzibb2.group, Mbb.group)

ggplot()+
  geom_point(data = viz[viz$gene_name %in% promising.genes, ],
             aes(x = gene_name, y = gene_usage_pct, fill = condition),
             shape = 21, stroke = 0, alpha = 0.5, size = 2,
             position = position_jitterdodge(jitter.width = 0.15,
                                             dodge.width = 0.85,
                                             jitter.height = 0))+
  geom_errorbar(data = group[group$gene_name %in% promising.genes, ],
                aes(x = gene_name, ymin = ppc.L,
                    ymax = ppc.H, linetype = condition, col = model),
                position = position_dodge(width = .8), width = 0.75)+
  theme_bw(base_size = 9)+
  theme(legend.position = "top")+
  scale_fill_manual(values = c("black", "darkgray"))










# 3-9
Mbb$usage.data$gene_names

e.bb <- rstan::extract(object = Mbb$glm)
e.zibb <- rstan::extract(object = Mzibb$glm)
e.zibb2 <- rstan::extract(object = Mzibb2$glm)
e.mix <- rstan::extract(object = Mmix$glm)

getGeneD <- function(e, gene_nr) {
  m1 <- e$alpha_gene[, gene_nr] + (1)*e$beta_gene[, gene_nr]
  p1 <- 1/(1 + exp(-(m1)))

  mn1 <- e$alpha_gene[, gene_nr] + (-1)*e$beta_gene[, gene_nr]
  pn1 <- 1/(1 + exp(-(mn1)))
  return (list(p1 = p1, pn1 = pn1))
}

e.post.bb <- getGeneD(e = e.bb, gene_nr = 2)
plot(density(e.post.bb$p1*100))
points(density(e.post.bb$pn1*100), col = "red", type = "l", xlim = c(0, 10))
plot(density(e.post.bb$pn1*100))


e.post.zibb <- getGeneD(e = e.zibb, gene_nr = 2)
plot(density(e.post.zibb$p1*100))
points(density(e.post.zibb$pn1*100), col = "red", type = "l", xlim = c(0, 10))
plot(density(e.post.zibb$pn1*100))



e.post.mix <- getGeneD(e = e.mix, gene_nr = 2)
plot(density(e.post.mix$p1*100))
points(density(e.post.mix$pn1*100), col = "red", type = "l", xlim = c(0, 10))
plot(density(e.post.mix$p1*100))


plot(density(e.zibb$Yhat_individual[1:1000, 2,1]))
