
mapq <- 30
pl <- "ILLUMINA"
myDefault <- list(mapq = mapq, pl = pl)
myArgs <- commandArgs(trailingOnly = T)
ind <- myArgs[1]
d <- myArgs[2]

path <- paste(getwd(), "/FASTQ/", ind, "/",sep="")
fastq <- list.files(path, pattern = ".fastq")

myData <- data.frame(ID = fastq, Data = fastq, MAPQ = mapq,
                    LB = fastq, PL = pl, SM = d)

myData$ID <- sub(paste(".*", ind, "_", sep = ""), "", myData$ID)
myData$ID <- sub(".fastq.gz", "", myData$ID)

myData$Data <- sub("^", path, myData$Data)

myData$LB <-  sub("_.*", "", sub(paste(".*",ind, "_", sep = ""), "", myData$LB))

write.table(myData, paste(d, ".txt", sep = ""), col.names = T,
            row.names = F, sep = "\t", quote = F)