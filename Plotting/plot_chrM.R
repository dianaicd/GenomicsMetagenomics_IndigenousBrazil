library(Gviz)
library(BSgenome.Hsapiens.UCSC.hg19)
sTrack <- SequenceTrack(Hsapiens)
afrom <- 0
ato <- 17000
alTrack <- AlignmentsTrack(range = "~/Projects/Botocudos/test/MA2392.hg19_mito.realigned.bam", 
                           isPaired = F)
plotTracks(sTrack, chromosome = "chrM", from = afrom, to = 1000)
png("~/Projects/Botocudos/Plots/mito_reads.png")
plotTracks(c(alTrack, sTrack), chromosome = "chrM", from = afrom,to = ato, cex = 3, min.height = 10)
dev.off()
