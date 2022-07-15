args <- commandArgs(TRUE)
# coordinates SuruiA_decoded.Summary.txt 
path_varSummary <- args[1]
# which ones have data: SuruiA_decoded.All_posterior_probs.txt
path_varAll <- args[2]

# state of interest: Australasian
state <- args[3]

prefix_out <- args[4]

varSummary <- read.table(path_varSummary, header = T, sep = "\t")
varAll <- read.table(path_varAll, header = T, sep = "\t")

chunksState <- varSummary[varSummary$state == state, ]

#block <- chunksState[1,]

chunk_to_bed <- function(block, varAll, n_block){
    chrom <- as.numeric(block$chrom[1])
    start <- as.numeric(block$start[1])
    end <- as.numeric(block$end[1])
    prob <- block$mean_prob


    variants_block <- varAll[varAll$chrom == chrom & varAll$start >= block$start & varAll$start < block$end,]

    obs_vars <- unlist(
        sapply(
            variants_block$variants[variants_block$variants != ""],
            function(v) strsplit(v, ",")
        )
    )

    small_bed <- data.frame(
        chr = chrom,
        pos = obs_vars,
        prob = prob,
        n_block = n_block
    )

    return(small_bed)
}

myBed <- do.call(
        rbind,
        lapply(
            1:nrow(chunksState),
            function(i)
            chunk_to_bed(block = chunksState[i,], varAll = varAll, n_block = i)
        )
    )


write.table(myBed, paste0(prefix_out, ".txt"), col.names=T, row.names=F, sep = "\t", quote = F)

myBed$start <- as.numeric(myBed$pos) - 1
myBed$end <- as.numeric(myBed$pos)
write.table(myBed[,c("chr", "start", "pos")],  paste0(prefix_out, ".bed"), col.names=F, row.names=F, sep = "\t", quote = F)