#!/software/R/3.4.2/bin/Rscript

# Script to merge SNPs from different individuals
# The input files must have 4 columns: chr, pos, ref, alt
#

library(plyr)
pattern <- "sites_.*txt"
files <- list.files(pattern = pattern)

sites <- data.frame()
N <- length(files)
j <- 0
for(f in files){
    j <- j +1
    print(paste("Individual", j, "out of", N))
    print(paste("Merging", f))
    toMerge <- read.table(f, stringsAsFactors = F)
    colnames(toMerge) <- c("chr", "pos", "ref", "alt")

    toMerge$pos <- paste(toMerge$chr, toMerge$pos, sep = "_")
    toMerge$chr <- NULL
    toMerge$variants <- paste(toMerge$ref, toMerge$alt, sep = "")

    if(dim(sites)[1] == 0){
        sites <- toMerge
        sites$variants <- paste0(sites$ref, sites$alt, sep = "")
        sites$ref <- NULL
        sites$alt <- NULL
        next
    }

    print("Sorting")
    sites <- sites[order(sites$pos),]
    toMerge <- toMerge[order(toMerge$pos),]

    print(paste("Number of sites in this individual:", dim(toMerge)[1]))

    overlapSites <-  which(toMerge$pos %in% sites$pos)
    overlapToMerge  <-  which(sites$pos %in% toMerge$pos)
    
    shared <- data.frame(variants = sites$variants[overlapToMerge], 
                    newRef = toMerge$ref[overlapSites],
                    newAlt = toMerge$alt[overlapSites])
    print(paste("Positions already present in the panel:", dim(shared)[1]))
    shared$variants <- as.character(shared$variants)
    toMerge$ref <- NULL
    toMerge$alt <- NULL
    # Are the new alt alleles already in the dataset?
    for(i in 2:3){
        newAlleles <- unlist(apply(as.matrix(shared), 1, 
                            function(x) grep(x[i], x[1]))) -1
        if(sum(newAlleles)){
            print(paste(sum(newAlleles), colnames(shared)[i], "alleles"))
            index <- which(newAlleles == 1)
            shared$variants[index] <- paste0(shared$variants[index],
                                    shared$newRef[index], paste = "")
        }
    }

    sites$variants[overlapToMerge] <- shared$variants
    # Reduce size of data frame to speed up merging
    toMerge <- toMerge[-overlapSites,] 
    sites <- join(sites, toMerge,  by = c("pos", "variants"), type = "full")
    print(paste("Current number of sites:", dim(sites)[1]))

}

name <- paste("merged_", dim(sites)[1], "sites_Africans.txt", sep = "")
sites$chr <- unlist(strsplit(sites$pos, split = "_"))[seq(1, 2*dim(sites)[1], 2)]
sites$pos <- unlist(strsplit(sites$pos, split = "_"))[seq(2, 2*dim(sites)[1], 2)]

sites$ref <- unlist(strsplit(sites$variants, split = ""))[seq(1, 2*dim(sites)[1], 2)]
sites$alt <- unlist(strsplit(sites$variants, split = ""))[seq(2, 2*dim(sites)[1], 2)]
sites$variants <- NULL

#head(sites[,c("chr", "pos", "ref", "alt")])

sites <- sites[,c("chr", "pos", "ref", "alt")]
sites <- sites[order(sites$chr, sites$pos),]
write.table(sites, name, col.names = F, 
            row.names = F, quote = F, sep = "\t")

