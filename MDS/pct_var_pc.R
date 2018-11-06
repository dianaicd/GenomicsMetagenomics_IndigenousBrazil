pct_var_pc <- function(eval, pc){
  pct.varPC = paste("PC",pc," (",round(100*eval[pc,]/sum(eval),2),"%)",sep="")
  return(pct.varPC)
}

