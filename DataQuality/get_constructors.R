
get_constructor <- function(my_table, r=1){
  msg <- "Error building the constructors"
  result <- tryCatch({ds.mincount.bootstrap(n = my_table, r = r)},
                     error =function(err){ warning(msg) }
  )
  if(result == msg){
    times <- 100
    while(result == msg & times >1){
      result <- tryCatch({
        ds.mincount.bootstrap(n = my_table, r = r, times = times)},
        error =function(err){ warning(msg) })
      times <- round(times/2)
    }
  }
  return(result)
}
