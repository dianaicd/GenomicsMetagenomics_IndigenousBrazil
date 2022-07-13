args <- commandArgs(TRUE)

## Default setting when no arguments passed
if(length(args) < 1) {
  args <- c("--help")
}

## Help section
if("--help" %in% args) {
  cat("
 
      Arguments:
    --input_sfs
    --pop1
    --pop2
    --swap
    --order_pop1
    --order_pop2
    --out
    --help                              - print this text
 
      Example:
Rscript rearange_2dsfs.R --input_sfs=pop1_pop2.sfs --pop1=pop1 --pop2=pop2 --swap=TRUE --order_pop1=1 --order_pop2=2 --out=new_pop2_pop1.sfs \n\n")
  
  q(save="no")
}

## Parse arguments (we expect the form --arg=value)
parseArgs <- function(all_args){
  lapply(
    all_args, 
    function(str) {
      # removing -- from the name; this will go to the first column of argsDF
      new_str = sub("^--", "", str)
      regmatches(new_str, regexpr("=", new_str), invert = TRUE)[[1]]
    }
    )
} 

argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1


get_args <- function(argsL, name, default=NA){
  if(name %in% names(argsL)){
    value = argsL[[name]]
  }else{
    value = default
  }
  return(value)
}


# getting arguments
filename          = get_args(argsL, "input_sfs", NA)
ppop1 = get_args(argsL, "pop1")
pop2 = get_args(argsL, "pop2")
swap = eval( get_args(argsL, "swap", FALSE) )
order_pop1 = get_args(argsL, "order_pop1", "0")
order_pop2 = get_args(argsL, "order_pop2", "1")
out = get_args(argsL, "out")

sfs <- read.table(filename)
new_sfs <- sfs

colnames(new_sfs) <- sub(
    "d._", paste0("d", order_pop1, "_"), colnames(new_sfs)
)
rownames(new_sfs) <- sub(
    "d._", paste0("d", order_pop2, "_"), rownames(new_sfs)
)
if(swap){
    new_sfs <- t(new_sfs)
}

write.table(new_sfs, out, col.names = T, row.names = T, sep = "\t", quote = F)