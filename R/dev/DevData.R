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
