require(ggplot2)
require(gridExtra)
require(ggrepel)
require(rstan)
require(tidyr)
require(reshape2)
rstan_options(auto_write = TRUE)

source(file = "R/Util.R")
source(file = "R/Usage.R")
data("IGHV_HCV")




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
viz$gene_name_fig <- gsub(pattern = "IGHV", replacement = '', x = viz$gene_name)
viz$condition <- ifelse(test = viz$condition == "hcv", yes = "HCV", no = "HD")


M <- get(load(file = "R/dev/ighv_hcv_zibb_model.RData"))




# format data
stats <- M$glm.summary
stats$effect_col <- ifelse(test = stats$effect_L <= 0 & stats$effect_H >= 0,
                           yes = "include", no = "exclude")
stats <- stats[order(abs(stats$effect_mean), decreasing = F), ]
stats$gene_fac <- factor(x = stats$gene_name,
                         levels = stats$gene_name)
stats <- merge(x= stats, y = M$test.stats, by = "gene_name")




stats$gene_name_fig <- gsub(pattern = "IGHV", replacement = '', x = stats$gene_name)
stats$t.test.fdr.pvalue[is.finite(stats$t.test.fdr.pvalue) == F] <- 1
stats$u.test.fdr.pvalue[is.finite(stats$u.test.fdr.pvalue) == F] <- 1




# T-test
g <- ggplot()+
  geom_hline(yintercept = c(-log10(0.05), -log10(0.01)), linetype = "dashed", col = "darkgray")+
  geom_text_repel(data = stats[stats$pmax >= 0.8, ], min.segment.length = 0.1, size = 2,
                  aes(x = pmax, y = -log10(t.test.fdr.pvalue), label = gene_name_fig))+
  geom_point(data = stats, aes(x = pmax, y = -log10(t.test.fdr.pvalue)),
             shape = 21, fill = "red", col = "black", size = 1.5)+
  xlim(0.5, 1)+
  ylab(label = "P-value [-log10]")+
  xlab(label = expression(pi))+
  theme_bw(base_size = 9)+
  theme(legend.position = "top")+
  scale_color_discrete(name = '')
ggsave(filename = "R/dev/manuscript/ttest.eps", plot = g,
       device = "eps", width = 4, height = 2.5, dpi = 600)







# U-test
g <- ggplot()+
  geom_hline(yintercept = c(-log10(0.05), -log10(0.01)), linetype = "dashed", col = "darkgray")+
  geom_text_repel(data = stats[stats$pmax >= 0.8, ], min.segment.length = 0.1, size = 2,
                  aes(x = pmax, y = -log10(u.test.fdr.pvalue), label = gene_name_fig))+
  geom_point(data = stats, aes(x = pmax, y = -log10(u.test.fdr.pvalue)),
             shape = 21, fill = "red", col = "black", size = 1.25)+
  xlim(0.5, 1)+
  ylab(label = "P-value [-log10]")+
  xlab(label = expression(pi))+
  theme_bw(base_size = 9)+
  theme(legend.position = "top")+
  scale_color_discrete(name = '')
ggsave(filename = "R/dev/manuscript/utest.eps", plot = g,
       device = "eps", width = 4, height = 2.5, dpi = 600)





group.ppc <- M$group.ppc$group.ppc
group.ppc <- group.ppc[group.ppc$gene_name %in% promising.genes, ]

point.data <- M$group.ppc$individual.pct.data
point.data <- point.data[point.data$gene_name %in% promising.genes, ]

ggplot()+
  geom_errorbar(data = group.ppc,
                aes(x = gene_name, ymin = ppc.L,
                    ymax = ppc.H, col = condition),
                position = position_dodge(width = .8), width = 0.75)+
  geom_point(data = point.data,
             aes(x = gene_name, y = gene_usage_pct, col = condition),
             shape = 21, size = 1.5, fill = "black",
             position = position_jitterdodge(jitter.width = 0.15,
                                             jitter.height = 0,
                                             dodge.width = 0.8))+
  theme_bw(base_size = 9)+
  theme(legend.position = "top")


promising.genes <- c("IGHV1-58", "IGHV3-72", "IGHV4-30-4", "IGHV3-30-3", "IGHV3-23", "IGHV3-7")


g <- ggplot(data = viz[viz$gene_name %in% promising.genes, ])+
  facet_wrap(facets = ~gene_name_fig, ncol = 2, scales = "free")+
  geom_point(aes(x = '', y = gene_usage_count, fill = condition, shape = condition, size = total/10^3),
             position = position_jitterdodge(dodge.width = 1, jitter.width = 0.35, jitter.height = 0),
             stroke = 0.5)+
  theme_bw(base_size = 9)+
  ylab(label = "Usage [count]")+
  xlab(label = '')+
  theme(legend.position = "top")+
  scale_fill_manual(name = "Condition", values = c("orange", "#a4c0ed"))+
  scale_shape_manual(name = "Condition", values = c(21, 22))+
  scale_size_continuous(name = "N (in thousands)", range = c(0.75, 2.5))+
  guides(size = guide_legend(nrow = 1, byrow = TRUE))


g
ggsave(filename = "R/dev/manuscript/results.eps", plot = g,
       device = "eps", width = 2.5, height = 5, dpi = 600)


