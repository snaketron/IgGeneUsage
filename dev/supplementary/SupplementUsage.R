require(ggplot2)
require(ggforce)


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


# viz$is.small <- viz$gene_usage_pct <= 2
# t <- table(viz$gene_name, viz$is.small)
# viz <- viz[!viz$gene_name %in% names(which(t[, 2] == 29)), ]
viz$gene_name_fig <- gsub(pattern = "IGHV", replacement = '', x = viz$gene_name)
viz$condition <- ifelse(test = viz$condition == "hcv", yes = "HCV", no = "HD")
viz$family <- substr(x = viz$gene_name_fig, start = 1, stop = 1)
viz$sample_nr <- as.numeric(as.factor(viz$sample_id))




g <- ggplot(data = viz)+
  facet_wrap(facets = ~gene_name, ncol = 7, scales = "free_y")+
  geom_point(aes(x = sample_nr, y = gene_usage_pct/100, fill = condition,
                 shape = condition), size = 1, stroke = 0.15)+
  theme_bw(base_size = 9)+
  ylab(label = "Usage")+
  xlab(label = '')+
  scale_fill_manual(name = "condition", values = c("orange", "#a4c0ed"))+
  scale_shape_manual(name = "condition", values = c(21, 22))+
  scale_x_continuous(breaks = c(1, 10, 20, 30), labels =  c(1, 10, 20, 30))+
  scale_y_continuous(labels = scales::scientific,
                     breaks = scales::pretty_breaks(n = 3))+
  # scale_y_continuous(breaks = scales::pretty_breaks(n = 3))+
  theme(legend.position = "none")
g

ggsave(filename = "dev/supplementary/UsageProportions.pdf", plot = g,
       device = "pdf", width = 7, height = 7.25, dpi = 600)
ggsave(filename = "dev/supplementary/UsageProportions.eps", plot = g,
       device = "eps", width = 7, height = 7.25, dpi = 600)





g <- ggplot(data = viz)+
  facet_wrap(facets = ~gene_name, ncol = 7, scales = "free_y")+
  geom_point(aes(x = sample_nr, y = gene_usage_count, fill = condition,
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

ggsave(filename = "dev/supplementary/UsageCount.pdf", plot = g,
       device = "pdf", width = 7, height = 7.25, dpi = 600)
ggsave(filename = "dev/supplementary/UsageCount.eps", plot = g,
       device = "eps", width = 7, height = 7.25, dpi = 600)
