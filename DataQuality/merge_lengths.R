merge_lengths <- function(summ, path){
  require(plyr)
  types <- c("", ".mito") # nuclear and mitochondrial
  individuals <- unique(boto$sample)
  for(ind in individuals){
    libs <- boto$library[boto$sample == ind & boto$library != "All" & boto$library != "*"]
    length_all <- data.frame()
    for(type in types){
      for(lib in libs){
        
        length_tmp <- read.table(paste(path, ind, "/", 
                                       lib, type, ".length", sep = ""))
        colnames(length_tmp) <- c("length", lib)
        if(dim(length_all)[1]){
          length_all <- join(length_all, length_tmp, by = "length")
        }else{
          length_all <- length_tmp
        }
        
      }
      length_all$All <- ifelse(ncol(length_all == 2), length_all[,2],
                               rowSums(length_all[,2:ncol(length_all)], na.rm = T))
      new_path <- paste(path,  ind, "/All", 
                        type, ".length", sep = "")
      print(paste("will save to ", new_path, sep = ""))
      write.table(length_all[,c(1, ncol(length_all))],
                  new_path, col.names = F, row.names = F, quote = F,
                  sep = "\t")
    }
    
  }
}
