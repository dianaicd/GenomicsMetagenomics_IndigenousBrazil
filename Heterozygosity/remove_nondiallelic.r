args <- commandArgs(TRUE)
nondiallelic <- args[1]
sample <- args[2]

toRemove <- read.table(nondiallelic)
all <- read.table(sample)

toRemove$sites <- paste(toRemove$V1, toRemove$V2, sep = "_")
all$sites <- paste(all$V1, all$V2, sep = "_")
nSNP <- dim(all)[1]

all <- all[!(all$sites %in% toRemove$sites),]
nRemoved <- nSNP - dim(all)[1]
all$sites <- NULL

write.table(all, sample, col.names = F, row.names = F, quote = F, sep = "\t")
print(paste("Removed", nRemoved, "sites from", sample))
