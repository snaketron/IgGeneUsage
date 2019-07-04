# Data preparation
IGHV_HCV <- get(load(file = "inst/BCR_HCV_usage.RData"))
IGHV_HCV <- IGHV_HCV[which(regexpr(pattern = "IGHV", text = IGHV_HCV$key) != -1), ]

# CSM
IGHV_HCV_CSM <- IGHV_HCV[IGHV_HCV$cell %in% "CSM" & IGHV_HCV$replicate == 1, ]
IGHV_HCV_CSM$sample_id <- as.character(IGHV_HCV_CSM$patient)
IGHV_HCV_CSM$condition <- IGHV_HCV_CSM$status
IGHV_HCV_CSM$gene_name <- IGHV_HCV_CSM$key
IGHV_HCV_CSM$gene_usage_count <- IGHV_HCV_CSM$Y
save(IGHV_HCV_CSM, file = "inst/IGHV_HCV_CSM.RData")
write.table(x = IGHV_HCV_CSM, file = "inst/IGHV_HCV_CSM.csv",
            row.names = F)

# naive
IGHV_HCV_naive <- IGHV_HCV[IGHV_HCV$cell %in% "naive", ]
IGHV_HCV_naive$sample_id <- as.character(IGHV_HCV_naive$patient)
IGHV_HCV_naive$condition <- IGHV_HCV_naive$status
IGHV_HCV_naive$gene_name <- IGHV_HCV_naive$key
IGHV_HCV_naive$gene_usage_count <- IGHV_HCV_naive$Y
save(IGHV_HCV_naive, file = "inst/IGHV_HCV_naive.RData")

# trbv ms
require(immunarch)
immunarch::gene_stats()
data(immdata)
imm_gu = geneUsage(immdata$data, "hs.trbv")

require(reshape2)
d <- melt(data = imm_gu, value.name = "gene_usage_count", na.rm = F)
colnames(d) <- c("gene_name", "sample_id", "gene_usage_count")
d$condition <- ifelse(test = regexpr(pattern = "MS", text = d$sample_id) != -1, yes = "MS", no = "Control")
d <- d[, c("sample_id", "condition", "gene_name", "gene_usage_count")]
d$gene_usage_count[is.na(d$gene_usage_count) == T] <- 0

write.table(x = d, file = "inst/TRBV_MS.csv",
            row.names = F)


# MS VDJ
# basic <- read.csv(file = "~/Desktop/work/tools/vdjtools-1.2.1/basic_u/basicstats.txt", sep = "\t", as.is = T)
# usage <- read.csv(file = "~/Desktop/work/tools/vdjtools-1.2.1/vusage_u/segments.unwt.V.txt", sep = "\t", as.is = T)
# t <- merge(x = basic, y = usage, by = c("sample_id", "sex", "age", "state", "lane"))
# is <- which(regexpr(pattern = "TRBV", text = colnames(t)) != -1)
# c <- t
# for(i in 1:length(is)) {
#   c[, is[i]] <- t[, is[i]]*t$diversity
# }
# c <- c[, c(1, is)]
# k <- melt(data = c, id.vars = "sample_id", value.name = "gene_usage_count")
# colnames(k) <- c("sample_id", "gene_name", "gene_usage_count")
# k$gene_name <- as.character(k$gene_name)
# k$gene_name <- gsub(pattern = '\\.', replacement = '-', x = k$gene_name)
# k$condition <- ifelse(test = regexpr(pattern = "C", text = k$sample_id) != -1, yes = "C", no = "MS")
# d <- k
# d <- d[, c("sample_id", "condition", "gene_name", "gene_usage_count")]
# write.table(x = d, file = "inst/TRBV_MS.csv", row.names = F)


# IgM at -1h, +7d
library(alakazam)
library(dplyr)
library(scales)
data(ExampleDb)
family.clone <-countGenes(ExampleDb,
                          gene="V_CALL",
                          groups=c("SAMPLE", "ISOTYPE"),
                          mode="family",
                          clone="CLONE")
d <- family.clone
d$sample_id <- paste(d$SAMPLE, d$ISOTYPE, sep = "_")
d$condition <- d$SAMPLE
d$gene_name <- d$GENE
d$gene_usage_count <- d$CLONE_COUNT
d <- d[, c("sample_id", "condition",
           "gene_name", "gene_usage_count")]
write.table(x = d, file = "inst/Ig.csv", row.names = F)

