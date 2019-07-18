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



g <- ggplot(data = viz)+
  facet_wrap(facets = ~gene_name_fig, ncol = 7, scales = "free")+
  geom_point(aes(x = gene_name_fig, y = gene_usage_pct, fill = condition,
                 shape = condition), size = 1, stroke = 0.15,
             position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.25))+
  theme_bw(base_size = 9)+
  ylab(label = "Usage [%]")+
  xlab(label = '')+
  theme(legend.position = "top")+
  scale_fill_manual(name = "condition", values = c("orange", "#a4c0ed"))+
  scale_shape_manual(name = "condition", values = c(21, 22))
g



ggsave(filename = "dev/supplementary/usage_pct.eps", plot = g,
       device = "eps", width = 6.5, height = 8, dpi = 600)





g <- ggplot(data = viz)+
  facet_wrap(facets = ~gene_name, ncol = 3, scales = "free")+
  geom_point(aes(x = sample_id, y = gene_usage_count, fill = condition,
                 shape = condition), size = 1, stroke = 0.15)+
  theme_bw(base_size = 9)+
  ylab(label = "Usage [count]")+
  xlab(label = '')+
  scale_fill_manual(name = "condition", values = c("orange", "#a4c0ed"))+
  scale_shape_manual(name = "condition", values = c(21, 22))+
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4))
g



ggsave(filename = "dev/supplementary/usage_count.eps", plot = g,
       device = "eps", width = 9, height = 24, dpi = 600)
