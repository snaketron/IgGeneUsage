individual <- read.csv(file = "~/Desktop/Tcell_data/individual__of_project_1561635663523.csv")
sample_descr <- read.csv(file = "~/Desktop/Tcell_data/sample_description__of_project_1561635690270.csv")
sample_profile <- read.csv(file = "~/Desktop/Tcell_data/sample_profile__of_project_1561635763679.csv")

usage.files <- list.files(path = "~/Desktop/Tcell_data/", full.names = T, pattern = "usage")
usage.files.short <- list.files(path = "~/Desktop/Tcell_data/", full.names = F, pattern = "usage")
usage.files.short <- do.call(rbind, strsplit(x = usage.files.short, split = "\\_jusage|\\_vusage"))[, 1]


usage.data <- c()
for(i in 1:length(usage.files)) {
  d <- read.csv(file = usage.files[i], as.is = T)
  d$SampleID <- usage.files.short[i]
  usage.data <- rbind(usage.data, d)
}
rm(d, i)
usage.data <- merge(x = usage.data, y = sample_profile, by = "SampleID")

save(usage.data, file = "T")
write.table(x = usage.data, file = "TCR_Cancer.csv", sep = ";", row.names = F)
