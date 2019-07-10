#### 1. General view ####
data(Ig)
head(Ig)

# we can compute the total number of rearrangements per sample
total.usage <- aggregate(gene_usage_count~sample_id+condition,
                         FUN = sum, data = Ig)
total.usage$total <- total.usage$gene_usage_count
total.usage$gene_usage_count <- NULL

# merge it with the original data
viz <- merge(x = Ig, y = total.usage,
             by = c("sample_id", "condition"),
             all.x = TRUE)

# compute %
viz$gene_usage_pct <- viz$gene_usage_count/viz$total*100

# visualize
ggplot(data = viz)+
  geom_point(aes(x = gene_name, y = gene_usage_pct, fill = condition,
                 shape = condition), stroke = 0, alpha = 0.5, size = 2,
             position = position_jitterdodge(jitter.width = 0.25,
                                             dodge.width = 0.85,
                                             jitter.height = 0))+
  theme_bw(base_size = 9)+
  ylab(label = "Usage [%]")+
  xlab(label = '')+
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4))+
  scale_fill_manual(name = "Condition", values = c("orange", "#6666ff"))+
  scale_shape_manual(name = "Condition", values = c(21, 22))



#### 2. Model parameters ####

# Mmix
# Mmix <- get(load(file = "R/dev/alakazam_ighv_families_zmix_model.RData"))
# Mmix.stats <- Mmix$glm.summary
# Mmix.stats$effect_col <- ifelse(test = Mmix.stats$effect_L <= 0 & Mmix.stats$effect_H >= 0,
#                                 yes = "include", no = "exclude")
# Mmix.stats <- Mmix.stats[order(abs(Mmix.stats$effect_mean), decreasing = F), ]
# Mmix.stats$gene_fac <- factor(x = Mmix.stats$gene_name,
#                               levels = Mmix.stats$gene_name)
# Mmix.stats <- merge(x = Mmix.stats, y = Mmix$test.stats, by = "gene_name")
# Mmix.stats$model <- "Mix"



# Mbin
Mbin <- get(load(file = "R/dev/alakazam_binomial.RData"))
Mbin.stats <- Mbin$glm.summary
Mbin.stats$effect_col <- ifelse(test = Mbin.stats$effect_L <= 0 & Mbin.stats$effect_H >= 0,
                               yes = "include", no = "exclude")
Mbin.stats <- Mbin.stats[order(abs(Mbin.stats$effect_mean), decreasing = F), ]
Mbin.stats$gene_fac <- factor(x = Mbin.stats$gene_name,
                             levels = Mbin.stats$gene_name)
Mbin.stats <- merge(x = Mbin.stats, y = Mbin$test.stats, by = "gene_name")
Mbin.stats$model <- "Binomial"


# Mbb
Mbb <- get(load(file = "R/dev/alakazam_beta_binomial_model.RData"))
Mbb.stats <- Mbb$glm.summary
Mbb.stats$effect_col <- ifelse(test = Mbb.stats$effect_L <= 0 & Mbb.stats$effect_H >= 0,
                                yes = "include", no = "exclude")
Mbb.stats <- Mbb.stats[order(abs(Mbb.stats$effect_mean), decreasing = F), ]
Mbb.stats$gene_fac <- factor(x = Mbb.stats$gene_name,
                              levels = Mbb.stats$gene_name)
Mbb.stats <- merge(x = Mbb.stats, y = Mbb$test.stats, by = "gene_name")
Mbb.stats$model <- "Beta-binomial"


# Mzibb
Mzibb <- get(load(file = "R/dev/alakazam_zibb_model.RData"))
Mzibb.stats <- Mzibb$glm.summary
Mzibb.stats$effect_col <- ifelse(test = Mzibb.stats$effect_L <= 0 & Mzibb.stats$effect_H >= 0,
                               yes = "include", no = "exclude")
Mzibb.stats <- Mzibb.stats[order(abs(Mzibb.stats$effect_mean), decreasing = F), ]
Mzibb.stats$gene_fac <- factor(x = Mzibb.stats$gene_name,
                             levels = Mzibb.stats$gene_name)
Mzibb.stats <- merge(x = Mzibb.stats, y = Mzibb$test.stats, by = "gene_name")
Mzibb.stats$model <- "ZIBB"


stats <- rbind(Mbin.stats, Mbb.stats, Mzibb.stats)
rm(Mbin.stats, Mbb.stats, Mzibb.stats)


# visualize disruption effects on each gene
ggplot(data = stats)+
  geom_vline(xintercept = 0, linetype = "dashed", col = "darkgray")+
  geom_errorbarh(aes(y = gene_fac, xmin = effect_L,
                     xmax = effect_H, col = model))+
  geom_point(aes(y = gene_fac, x = effect_mean, col = model))+
  theme_bw(base_size = 9)+
  theme(legend.position = "top")+
  xlab(label = "Mean effect")






####  LOO ####

loo::loo(Mbin$glm)

loo::loo(Mbb$glm)

loo::loo(Mzibb$glm)


loo::compare(loo::loo(Mbin$glm),
             loo::loo(Mbb$glm),
             loo::loo(Mzibb$glm))




#### 3. Model predictions (gene level) ####
# Mmix.individual <- Mmix$ppc
# Mmix.individual$model <- "Mix"

Mbin.individual <- Mbin$ppc
Mbin.individual$model <- "Binomial"

Mzibb.individual <- Mzibb$ppc
Mzibb.individual$model <- "ZIBB"

Mbb.individual <- Mbb$ppc
Mbb.individual$model <- "Beta-binomial"

individual <- rbind(Mbin.individual, Mbb.individual, Mzibb.individual)
rm(Mbin.individual, Mbb.individual, Mzibb.individual)



ggplot()+
  facet_grid(facets = model~gene_name, scales = "free_x")+
  geom_errorbar(data = individual,
                aes(x = sample_id, ymin = ppc.pct.L,
                    ymax = ppc.pct.H, col = condition),
                position = position_dodge(width = .8), width = 0.75)+
  geom_point(data = viz,
             aes(x = sample_id, y = gene_usage_pct, fill = condition,
                 shape = condition))+
  theme_bw(base_size = 9)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4))+
  theme(legend.position = "top")






#### 4. Model predictions (individual level) ####
# Mmix.group <- Mmix$group.ppc$group.ppc
# Mmix.group$model <- "Mix"

Mbin.group <- Mbin$group.ppc$group.ppc
Mbin.group$model <- "Binomial"

Mzibb.group <- Mzibb$group.ppc$group.ppc
Mzibb.group$model <- "ZIBB"

Mbb.group <- Mbb$group.ppc$group.ppc
Mbb.group$model <- "Beta-binomial"

group <- rbind(Mbin.group, Mzibb.group, Mbb.group)
rm(Mbin.group, Mzibb.group, Mbb.group)



ggplot()+
  geom_point(data = viz, aes(x = gene_name, y = gene_usage_pct, fill = condition),
             shape = 21, stroke = 0, alpha = 0.5, size = 2,
             position = position_jitterdodge(jitter.width = 0.15,
                                             dodge.width = 0.85,
                                             jitter.height = 0))+
  geom_errorbar(data = group, aes(x = gene_name, ymin = ppc.L,
                    ymax = ppc.H, linetype = condition, col = model),
                position = position_dodge(width = .8), width = 0.75)+
  theme_bw(base_size = 9)+
  theme(legend.position = "top")+
  scale_fill_manual(values = c("black", "darkgray"))


