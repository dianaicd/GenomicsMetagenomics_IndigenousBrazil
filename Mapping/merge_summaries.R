setwd("/scratch/beegfs/monthly/dcruzdva/Botocudos/BAM/2018_10_23/test")
library(plyr)
files <- list.files()
index <- grepl(".summary$", files)
files <- files[index]
boto <- data.frame()
for(f in files){
  tmp <- read.table(f, header = T)
  ind <- gsub("\\..*", "", f)
#  histo <- read.table(paste("./", ind, "_histogram.coverage", sep = ""))
#  bases <- sum(histo$V3[histo$V1 == "genome" & histo$V2 != 0])
#  tmp$Bases_covered <- bases
#  tmp$Bases_covered_fraction <- bases/sum(histo$V3[histo$V1 == "genome"])
  boto <- rbind(boto, tmp)
}

write.table(boto,
            paste("Botocudos_summary_", Sys.Date(), ".table", sep = ""),
            row.names = F, col.names = T, quote = F)
