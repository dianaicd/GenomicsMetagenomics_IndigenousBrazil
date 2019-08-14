require(preseqR)
require(parallel)

# Constructors
# A model is fitted to our data. We will have a constructor 
# that we can use later to estimate throughput for
# larger experiments
get_constructor <- function(my_table, r=1){
  msg <- "Error building the constructors"
  result <- tryCatch({ds.rSAC.bootstrap(n = my_table, r = r)},
                     error =function(err){ warning(msg) }
  )
  if(result == msg){
    times <- 100
    while(result == msg){
      result <- tryCatch({
        ds.rSAC.bootstrap(n = my_table, r = r, times = times)},
        error =function(err){ warning(msg) })
      #times <- round(times/2)
    }
  }
  return(result)
}
