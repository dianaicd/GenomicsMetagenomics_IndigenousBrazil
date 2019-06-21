library(plyr)
args = commandArgs(trailingOnly=TRUE)

panelName <- args[1]
bamlistName <- args[2]

panel <- read.table(paste(panelName, ".sites", sep = ""), header = F)
colnames(panel) <- c("chr", "pos", "major", "minor")
bamlist <- read.table(bamlistName, header = F)
individuals <- gsub(".*/", "", bamlist$V1)
individuals <- gsub(".bam", "", individuals)

x <- panel
for(ind in individuals){
    print(ind)
    y <- read.table(paste(ind, "/", ind, ".haplo.gz", sep = ""), header = T)
    y <- data.frame(chr = y$chr, pos = y$pos, ind = y$ind0)
    colnames(y) <- c("chr", "pos", ind)
    x <- join(x, y, by = c("chr", "pos"))
}

x[is.na(x),] <- 0