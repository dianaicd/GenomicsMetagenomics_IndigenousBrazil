options(scipen=999)
args = commandArgs(trailingOnly=TRUE)
blockSize <- as.integer(args[1])

x <- read.table("chrom_size.txt")
colnames(x) <- c("chr", "size")

blockStart <- 0    

for(chr in x$chr){
    ChrSize <- x$size[chr]

    while(blockStart < ChrSize){
        
        
        if(blockStart + blockSize < ChrSize){
            region <- data.frame("chr"=chr, "start"=blockStart,
                                "end"=blockStart + blockSize)
        
            fname <- paste("region_", chr, ":", blockStart, 
                            "-", blockStart + blockSize, ".txt", sep = "")
        }else{
           region <- data.frame("chr"=chr, "start"=blockStart,
                                "end"=ChrSize)
            if(chr+1 %in% x$chr){
                end <- signif(blockSize - ChrSize + blockStart, 1)
                newChr <- data.frame("chr"=chr+1, "start"=0,
                                "end"=end)
                region <- rbind(region, newChr)
                fname <- paste("region_", chr, ":", blockStart, 
                            "-", ChrSize,
                            "_", chr+1, ":", "0", "-",end ,".txt", sep = "")
            }else{
                fname <- paste("region_", chr, ":", blockStart, 
                            "-", ChrSize, ".txt", sep = "")
            }   
        }
        write.table(region, fname, col.names = F, row.names = F, quote = F,
                    sep = "\t")
        blockStart <- blockStart +blockSize
    }
    blockStart <- end
}