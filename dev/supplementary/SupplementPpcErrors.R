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




M <- get(load(file = "dev/ighv_hcv_zibb_model.RData"))




#### Error Repertoire: % ####
ppc <- M$ppc.data$ppc.repertoire
ppc$gene_name_fig <- gsub(pattern = "IGHV", replacement = '', x = ppc$gene_name)
ppc$condition <- ifelse(test = ppc$condition == "hcv", yes = "HCV", no = "HD")
ppc$family <- substr(x = ppc$gene_name_fig, start = 1, stop = 1)

g <- ggplot(data = ppc)+
  geom_sina(aes(x = condition, y = error.mean, fill = condition), shape = 21)+
  theme_bw(base_size = 9)+
  ylab(label = "Error [%]")+
  xlab(label = "Condition")+
  scale_fill_manual(name = "condition", values = c("orange", "#a4c0ed"))+
  theme(legend.position = "none")
g


ggsave(filename = "dev/supplementary/PpcGroupError.pdf", plot = g,
       device = "pdf", width = 4, height = 3, dpi = 600)
ggsave(filename = "dev/supplementary/PpcGroupError.eps", plot = g,
       device = "eps", width = 4, height = 3, dpi = 600)





#### Error Group: % ####
ppc <- M$ppc.data$ppc.gene
ppc$gene_name_fig <- gsub(pattern = "IGHV", replacement = '', x = ppc$gene_name)
ppc$condition <- ifelse(test = ppc$condition == "hcv", yes = "HCV", no = "HD")
ppc$family <- substr(x = ppc$gene_name_fig, start = 1, stop = 1)

g <- ggplot(data = ppc)+
  geom_sina(aes(x = condition, y = error.mean, fill = condition), shape = 21)+
  theme_bw(base_size = 9)+
  ylab(label = "Error [%]")+
  xlab(label = "Condition")+
  scale_fill_manual(name = "condition", values = c("orange", "#a4c0ed"))+
  theme(legend.position = "none")
g


ggsave(filename = "dev/supplementary/PpcGroupError.pdf", plot = g,
       device = "pdf", width = 4, height = 3, dpi = 600)
ggsave(filename = "dev/supplementary/PpcGroupError.eps", plot = g,
       device = "eps", width = 4, height = 3, dpi = 600)
