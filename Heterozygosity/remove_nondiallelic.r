args <- commandArgs(TRUE)
nondiallelic <- args[1]
sample <- args[2]

toRemove <- read.table(nondiallelic)
all <- read.table(sample)

toRemove$sites <- paste(toRemove$V1, toRemove$V2, sep = "_")
all$sites <- paste(all$V1, all$V2, sep = "_")
nSNP <- dim(all)[1]

all <- all[!(all$sites %in% toRemove$sites),]
nRemoved <- dim(all)[1]

print(paste("Removed", nRemoved, "sites from" sample))
