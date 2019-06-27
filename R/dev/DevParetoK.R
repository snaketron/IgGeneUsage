M3 <- get(load(file = "R/dev/ighv_hcv_beta_binomial_model.RData"))
rm(M)




loo.M3 <- loo::loo(x = loo::extract_log_lik(stanfit = M3$glm))

bad.k <- which(loo.M3$diagnostics$pareto_k >= 0.7)
dim(M3$usage.data$Y)

g <- floor(bad.k / 69)+1
q <- data.frame(sample.nr = g, bad.i = bad.k)
q$gene.rn <- q$bad.i-(q$sample.nr-1)*69


q$Y <- NA
q$N <- NA
q$p <- NA
q$gene_name <- NA
for(i in 1:nrow(q)) {
  q$Y[i] <- M3$usage.data$Y[q$gene.rn[i], q$sample.nr[i]]
  q$N[i] <- M3$usage.data$N[q$sample.nr[i]]
  q$p[i] <- q$Y[i]/q$N[i]*100
  q$gene_name[i] <- M3$usage.data$gene_names[q$gene.rn[i]]
}
rm(i, g)



gene_name <- M3$usage.data$gene_names[as.numeric(names(
  sort(table(q$gene.rn), decreasing = T)[1:10]))]



# p <- bad.k[which(bad.k %% 69 == 26)]
# M <- M3
# q$
gene_name


sort(table(q$gene.rn), decreasing = T)[1:10]





group.ppc <- M$group.ppc$group.ppc
group.ppc <- group.ppc[group.ppc$gene_name %in% gene_name, ]
group.ppc$gene_name <- factor(group.ppc$gene_name, levels = gene_name)

point.data <- M$group.ppc$individual.pct.data
point.data <- point.data[point.data$gene_name %in% gene_name, ]
point.data$gene_name <- factor(point.data$gene_name, levels = gene_name)

q <- q[q$gene_name %in% gene_name, ]

ggplot()+
  geom_errorbar(data = group.ppc,
                aes(x = gene_name, ymin = ppc.L,
                    ymax = ppc.H, col = condition),
                position = position_dodge(width = .8), width = 0.75)+
  geom_point(data = point.data,
             aes(x = gene_name, y = gene_usage_pct, col = condition),
             shape = 21, size = 1.5, fill = "black",
             position = position_jitterdodge(jitter.width = 0.25,
                                             jitter.height = 0,
                                             dodge.width = 0.8))+
  geom_point(data = q, aes(x = gene_name, y = p), col = "blue", size = 1.5)+
  theme_bw(base_size = 9)+
  theme(legend.position = "top")


