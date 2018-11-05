col2pop <- function(colors, pop, n, k = 10){
  
  rownames(pop) <- seq(1,n)
  #indexes to return
  indexes <- c()
  for(p in unique(pop$population)){
    # print(p)
    # name the individuals in pop p
    index <- as.integer(rownames(pop[pop$population==p,])); 
    #print(c("index:",index))
    # maximum fraction of admixture in pop p
    x <- max(admix[, index]) ;
    #print(c("x:",x))
    # Are there more than 1 ind. in pop p?
    if(length(index)>1){
      # Among the K ancestral populations,
      # find the maximum ancestry prop. for pop p
      maxes <- apply(admix[,index], MARGIN = 1 ,  max); 
      #print(c("maxes:",maxes))
      i <- 1
      # The first color corresponds to the ancestry
      # found in the biggest proportion on pop p
      col_i <- which(t(maxes) == x); 
      #print(c("col_i",col_i))
      # is the index of this color already in the
     #indexes to return?  
      while(col_i[i] %in% indexes && i < length(col_i)){
        i <- i + 1; 
        #print(c("i",i))
      }
      if(!col_i[i] %in% indexes){
        indexes <- c(indexes, col_i[i]); 
        #print(c("indexes",indexes))
      }
      
    }else{
      maxes <- max(admix[,index]); 
      #maxes
      col_i <- which(t(maxes) == x); 
      #col_i
      if(!col_i %in% indexes){
        indexes <- c(indexes, col_i); 
        #indexes
      }
    }
  }
  left <- !seq(1,k) %in% indexes
  indexes <- c(indexes, seq(1, k)[left])
  return(indexes)
}



