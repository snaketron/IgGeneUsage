require(IgGeneUsage)
require(knitr)
require(ggplot2)
require(ggforce)
require(gridExtra)
require(ggrepel)
require(rstan)
require(reshape2)
rstan_options(auto_write = TRUE)




data("IGHV_HCV", package = "IgGeneUsage")



## Data visualization
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

# visualize
ggplot(data = viz)+
  geom_point(aes(x = gene_name, y = gene_usage_pct,
                 fill = condition, shape = condition),
             position = position_dodge(width = .7), stroke = 0)+
  theme_bw(base_size = 12)+
  ylab(label = "Usage [%]")+
  xlab(label = '')+
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4))+
  scale_fill_manual(name = "Condition", values = c("orange", "#6666ff"))+
  scale_shape_manual(name = "Condition", values = c(21, 22))


M <- get(load(file = "dev/ighv_hcv_zibb_model.RData"))




#### PPC: counts ####
ppc <- M$ppc.data$ppc.repertoire
ppc$gene_name_fig <- gsub(pattern = "IGHV", replacement = '', x = ppc$gene_name)
ppc$condition <- ifelse(test = ppc$condition == "hcv", yes = "HCV", no = "HD")
ppc$family <- substr(x = ppc$gene_name_fig, start = 1, stop = 1)
ppc$sample_nr <- as.numeric(as.factor(ppc$sample_id))

g <- ggplot(data = ppc)+
  facet_wrap(facets = ~gene_name, ncol = 7, scales = "free_y")+
  geom_errorbar(aes(x = sample_nr, ymin = ppc.count.L, ymax = ppc.count.H),
                col = "darkgray", width = 0)+
  geom_point(aes(x = sample_nr, y = observed.count, fill = condition,
                 shape = condition), size = 1, stroke = 0.15)+
  theme_bw(base_size = 9)+
  ylab(label = "Usage [count]")+
  xlab(label = '')+
  scale_fill_manual(name = "condition", values = c("orange", "#a4c0ed"))+
  scale_shape_manual(name = "condition", values = c(21, 22))+
  scale_x_continuous(breaks = c(1, 10, 20, 30), labels =  c(1, 10, 20, 30))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 3))+
  theme(legend.position = "none")
g

ggsave(filename = "dev/supplementary/PpcCount.pdf", plot = g,
       device = "pdf", width = 7, height = 7.25, dpi = 600)
ggsave(filename = "dev/supplementary/PpcCount.eps", plot = g,
       device = "eps", width = 7, height = 7.25, dpi = 600)




#### PPC: % ####
ppc <- M$ppc.data$ppc.repertoire
ppc$gene_name_fig <- gsub(pattern = "IGHV", replacement = '', x = ppc$gene_name)
ppc$condition <- ifelse(test = ppc$condition == "hcv", yes = "HCV", no = "HD")
ppc$family <- substr(x = ppc$gene_name_fig, start = 1, stop = 1)
ppc$sample_nr <- as.numeric(as.factor(ppc$sample_id))

g <- ggplot(data = ppc)+
  facet_wrap(facets = ~gene_name, ncol = 7, scales = "free_y")+
  geom_errorbar(aes(x = sample_nr, ymin = ppc.pct.L, ymax = ppc.pct.H),
                col = "darkgray", width = 0)+
  geom_point(aes(x = sample_nr, y = observed.pct, fill = condition,
                 shape = condition), size = 1, stroke = 0.15)+
  theme_bw(base_size = 9)+
  ylab(label = "Usage [%]")+
  xlab(label = '')+
  scale_fill_manual(name = "condition", values = c("orange", "#a4c0ed"))+
  scale_shape_manual(name = "condition", values = c(21, 22))+
  scale_x_continuous(breaks = c(1, 10, 20, 30), labels =  c(1, 10, 20, 30))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 3))+
  theme(legend.position = "none")
g

ggsave(filename = "dev/supplementary/PpcProp.pdf", plot = g,
       device = "pdf", width = 7, height = 7.25, dpi = 600)
ggsave(filename = "dev/supplementary/PpcProp.eps", plot = g,
       device = "eps", width = 7, height = 7.25, dpi = 600)




#### PPC Group: % ####
ppc <- M$ppc.data$ppc.gene
ppc$gene_name_fig <- gsub(pattern = "IGHV", replacement = '', x = ppc$gene_name)
ppc$condition <- ifelse(test = ppc$condition == "hcv", yes = "HCV", no = "HD")
ppc$family <- substr(x = ppc$gene_name_fig, start = 1, stop = 1)

g <- ggplot(data = ppc)+
  geom_abline(slope = 1, intercept = 0)+
  geom_errorbar(aes(x = observed.mean, ymin = ppc.L, ymax = ppc.H),
                col = "darkgray", width = 0)+
  geom_point(aes(x = observed.mean, y = ppc.mean, fill = condition,
                 shape = condition), size = 1, stroke = 0.1)+
  theme_bw(base_size = 9)+
  ylab(label = "Mean PPC [%]")+
  xlab(label = "Mean observed [%]")+
  scale_fill_manual(name = "condition", values = c("orange", "#a4c0ed"))+
  scale_shape_manual(name = "condition", values = c(21, 22))+
  theme(legend.position = "top")
g


ggsave(filename = "dev/supplementary/PpcGroup.pdf", plot = g,
       device = "pdf", width = 4, height = 4.25, dpi = 600)
ggsave(filename = "dev/supplementary/PpcGroup.eps", plot = g,
       device = "eps", width = 4, height = 4.25, dpi = 600)

