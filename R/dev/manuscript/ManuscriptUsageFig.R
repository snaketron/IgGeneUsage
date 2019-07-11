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


viz$is.small <- viz$gene_usage_pct <= 2
t <- table(viz$gene_name, viz$is.small)
viz <- viz[!viz$gene_name %in% names(which(t[, 2] == 29)), ]
viz$gene_name_fig <- gsub(pattern = "IGHV", replacement = '', x = viz$gene_name)
viz$condition <- ifelse(test = viz$condition == "hcv", yes = "HCV", no = "HD")



# OLD
# g <- ggplot(data = viz)+
#   geom_point(aes(x = gene_name_fig, y = gene_usage_pct, fill = condition, shape = condition),
#              size = 1, stroke = 0, position = position_jitterdodge(dodge.width = 0.8,
#                                                                    jitter.width = 0.25,
#                                                                    jitter.height = 0))+
#   theme_bw(base_size = 9)+
#   ylab(label = "Usage [%]")+
#   xlab(label = '')+
#   theme(legend.position = "top",
#         axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
#   scale_fill_manual(name = "Condition", values = c("orange", "#6666ff"))+
#   scale_shape_manual(name = "Condition", values = c(21, 22))




g <- ggplot(data = viz)+
  geom_sina(aes(x = gene_name_fig, y = gene_usage_pct, fill = condition, shape = condition),
             size = 1, stroke = 0.25, position = position_dodge(width = 0.8), adjust = 15)+
  theme_bw(base_size = 9)+
  ylab(label = "Usage [%]")+
  xlab(label = '')+
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  scale_fill_manual(name = "Condition", values = c("orange", "#a4c0ed"))+
  scale_shape_manual(name = "Condition", values = c(21, 22))


g
ggsave(filename = "R/dev/manuscript/usage.eps", plot = g,
       device = "eps", width = 7, height = 2.5,
       dpi = 600)
