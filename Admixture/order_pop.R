# Order populations in ancestry plot

order_pop <- function(pops = c("Botocudo", "SouthAmerica", "NorthAmerica", "Oceania"),
                       inds){
  index <- c()
  
  s <- seq(1, length(inds$Region))
  for(p in pops){
    
    tmp <- s[inds$Region %in% p]
    index <- c(index, tmp)
  }
  return(index)
}
