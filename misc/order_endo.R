order_endo <- function(Min = 0, Max = 0.8, boto){
  ordered <- data.frame(Library = boto$Library[
    order(boto$hits_unique_frac_endogenous)],
    Endogenous = boto$hits_unique_frac_endogenous[
      order(boto$hits_unique_frac_endogenous)]
  )
  ordered <- ordered[ordered$Library != "*",]
  ordered <- ordered[ordered$Endogenous >= Min & ordered$Endogenous <= Max,]
  return(ordered)
}