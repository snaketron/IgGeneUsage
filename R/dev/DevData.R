# Data preparation
IGHV_HCV <- get(load(file = "inst/BCR_HCV_usage.RData"))
IGHV_HCV <- IGHV_HCV[which(regexpr(pattern = "IGHV", text = IGHV_HCV$key) != -1), ]
IGHV_HCV <- IGHV_HCV[IGHV_HCV$cell %in% c("naive", "CSM"), ]

IGHV_HCV$sample_id <- as.character(IGHV_HCV$patient)
IGHV_HCV$condition <- IGHV_HCV$status
IGHV_HCV$gene_name <- IGHV_HCV$key
IGHV_HCV$gene_usage_count <- IGHV_HCV$Y
save(IGHV_HCV, file = "inst/IGHV_HCV.RData")
