library(ggplot2)
library(cowplot)
source("~/Projects/Botocudos/Scripts/translate_ids.R")

# ChrM
setwd("~/Projects/Botocudos/Files/Contamination/2017_09_15/Mito/")

alleles <- read.table(gzfile("counts_chrM.counts.gz"), header = T)
pos <- read.table(gzfile("counts_chrM.pos.gz"), header = T)
individuals <- read.table("bam_mito.filelist")
individuals$ma <- gsub("\\..*", "", individuals$V1)
individuals$mn <- sapply(individuals$ma, function(x) ma2mn(x))
individuals$Contamination <- 0
for(ind in seq(1,22)){
  total <- rowSums(alleles[((ind-1)*4+1):(ind*4)])
  consensus <- apply(alleles[((ind-1)*4+1):(ind*4)], 1, function(x) max(x))
  contaminants <- total - consensus
  cont <- sum(contaminants)/sum(total)
  individuals$Contamination[ind] <- cont
  individuals$positions[ind] <- sum(total != 0)
  individuals$bases[ind] <- sum(total)
  #print(cont)
}

ggplot(individuals, aes(x = mn, y = Contamination)) +
  geom_bar(stat = "identity", fill = "orange") +
  labs(x = NULL, y = "Contamination (fraction)", 
       title = "Contamination from mtDNA") +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_text(
    aes(x = mn, y = 0.01, 
        label = round(positions/16569, 2))) +
  geom_text(
    aes(x = mn, y = 0.012, 
        label = round(bases/positions, 2)))

# ChrX
# setwd("~/Projects/Botocudos/Files/Contamination/2017_09_15/Nuc/")
# alleles <- read.table(gzfile("counts_chrX.counts.gz"), header = T)
# pos <- read.table(gzfile("counts_chrX.pos.gz"), header = T)
# individuals <- read.table("bam_nuc.filelist")
# for(ind in seq(0,21)){
#   total <- rowSums(alleles[(ind*4+1):(ind*4+4)])
#   consensus <- apply(alleles[(ind*4+1):(ind*4+4)], 1, function(x) max(x))
#   contaminants <- total - consensus
#   
#   print(sum(contaminants)/sum(total))
# }

