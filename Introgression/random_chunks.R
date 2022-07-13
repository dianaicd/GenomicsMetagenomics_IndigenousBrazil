path <- "~/Projects/Botocudos/Files/findIntrogression/2022_03/"

x <- read.table("~/Projects/Botocudos/Files/findIntrogression/2022_03/SuruiA_decoded.Summary.txt", header = T)

chunks <- x[x$state == "Australasian",]

chr_lengths <- sapply(1:22, function(i) tail(x$end[x$chrom == i],1))
genome_len <- sum(chr_lengths)


sample_in_chr <- function(chr, chunk_length, chr_lengths){
  overlaps_introgressed <- T
  
  while(overlaps_introgressed){
    start <- sample(1:(chr_lengths[chr] - chunk_length), 1)
    end <- start + chunk_length
    
    if(
    # start of random chunk is included in introgressed segment
      sum(start > chunks$start[chunks$chrom == chr] & start < chunks$end[chunks$chrom == chr]) |
    # end of random chunk is included in introgressed segment
      sum(end > chunks$start[chunks$chrom == chr] & end < chunks$end[chunks$chrom == chr]) |
    # introgressed chunk is in random chunk
      sum(chunks$start[chunks$chrom == chr] >= start & chunks$end[chunks$chrom == chr] <= end)
      ){
        overlaps_introgressed <- T
    }else{
        overlaps_introgressed <- F
    }
    
    return(start)
  }
}

random_lengths <-  sample(chunks$length, nrow(chunks))

start_pos <- sapply(
  1:nrow(chunks),
  function(i)
    sample_in_chr(chr = chunks$chrom[i],
                  chunk_length = random_lengths[i],
                  chr_lengths = chr_lengths)
)



random_chunks <- data.frame(
  chr = chunks$chrom,
  start_pos = start_pos,
  end_pos = start_pos + random_lengths
)


random_chunks <- random_chunks[order(random_chunks$chr, random_chunks$start_pos),]
hist(random_chunks$end_pos - random_chunks$start_pos, breaks = 100)


write.table(
  random_chunks,
  paste0(path, "random_chunks.bed"),
  col.names = F,
  row.names = F,
  sep = "\t"
            )

