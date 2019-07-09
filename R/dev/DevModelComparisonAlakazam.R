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
Mmix <- get(load(file = "R/dev/alakazam_ighv_families_zmix_model.RData"))
Mmix.stats <- Mmix$glm.summary
Mmix.stats$effect_col <- ifelse(test = Mmix.stats$effect_L <= 0 & Mmix.stats$effect_H >= 0,
                                yes = "include", no = "exclude")
Mmix.stats <- Mmix.stats[order(abs(Mmix.stats$effect_mean), decreasing = F), ]
Mmix.stats$gene_fac <- factor(x = Mmix.stats$gene_name,
                              levels = Mmix.stats$gene_name)
Mmix.stats <- merge(x = Mmix.stats, y = Mmix$test.stats, by = "gene_name")
Mmix.stats$model <- "Mix"


# Mbb
Mbb <- get(load(file = "R/dev/alakazam_ighv_families_beta_binomial_model.RData"))
Mbb.stats <- Mbb$glm.summary
Mbb.stats$effect_col <- ifelse(test = Mbb.stats$effect_L <= 0 & Mbb.stats$effect_H >= 0,
                                yes = "include", no = "exclude")
Mbb.stats <- Mbb.stats[order(abs(Mbb.stats$effect_mean), decreasing = F), ]
Mbb.stats$gene_fac <- factor(x = Mbb.stats$gene_name,
                              levels = Mbb.stats$gene_name)
Mbb.stats <- merge(x = Mbb.stats, y = Mbb$test.stats, by = "gene_name")
Mbb.stats$model <- "Beta-binomial"


stats <- rbind(Mmix.stats, Mbb.stats)


# visualize disruption effects on each gene
ggplot(data = stats)+
  geom_vline(xintercept = 0, linetype = "dashed", col = "darkgray")+
  geom_errorbarh(aes(y = gene_fac, xmin = effect_L,
                     xmax = effect_H, col = model))+
  geom_point(aes(y = gene_fac, x = effect_mean, col = model))+
  theme_bw(base_size = 9)+
  theme(legend.position = "top")+
  xlab(label = "Mean effect")


# format data
stats <- Mbb$glm.summary
stats$effect_col <- ifelse(test = stats$effect_L <= 0 & stats$effect_H >= 0,
                           yes = "include", no = "exclude")
stats <- stats[order(abs(stats$effect_mean), decreasing = F), ]
stats$gene_fac <- factor(x = stats$gene_name,
                         levels = stats$gene_name)

stats <- merge(x = stats, y = M$test.stats, by = "gene_name")

# visualize disruption effects on each gene
ggplot(data = stats)+
  geom_vline(xintercept = 0, linetype = "dashed", col = "darkgray")+
  geom_errorbarh(aes(y = gene_fac, xmin = effect_L,
                     xmax = effect_H, col = effect_col))+
  geom_point(aes(y = gene_fac, x = effect_mean))+
  theme_bw(base_size = 9)+
  theme(legend.position = "top")+
  scale_color_manual(name = '', values = c("exclude"="red",
                                           "include"="black"))+
  xlab(label = "Mean effect")




#### 3. Model predictions (gene level) ####
Mmix.group <- Mmix$group.ppc$group.ppc
Mmix.group$model <- "Mix"
Mbb.group <- Mbb$group.ppc$group.ppc
Mbb.group$model <- "Beta-binomial"
group <- rbind(Mmix.group, Mbb.group)


ggplot()+
  geom_errorbar(data = group,
                aes(x = gene_name, ymin = ppc.L,
                    ymax = ppc.H, linetype = condition, col = model),
                position = position_dodge(width = .8), width = 0.75)+
  geom_point(data = viz,
             aes(x = gene_name, y = gene_usage_pct, fill = condition),
             shape = 21, stroke = 0, alpha = 0.5, size = 2,
             position = position_jitterdodge(jitter.width = 0.25,
                                             dodge.width = 0.85,
                                             jitter.height = 0))+
  theme_bw(base_size = 9)+
  theme(legend.position = "top")






#### 4. Model predictions (individual level) ####
Mmix.individual <- Mmix$ppc
Mmix.individual$model <- "Mix"
Mbb.individual <- Mbb$ppc
Mbb.individual$model <- "Beta-binomial"
individual <- rbind(Mmix.individual, Mbb.individual)


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





loo::loo(Mmix$glm)
loo::loo(Mbb$glm)
