# Code to parse argument in a code, originally written by me for the script to infer sex in mapache
require(stringr)
#args <- commandArgs(TRUE)

## Default setting when no arguments passed
if(length(args) < 1) {
  args <- c("--help")
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


get_args <- function(name, default = NA, mandatory = FALSE, eval_string = FALSE, verbose = FALSE){
  if(name %in% names(argsL)){
    value = argsL[[name]]
  }else if(mandatory){
    error_message = str_glue("Error: parameter {name} must be specified to run this script.")
    cat(error_message)
    quit(save = "no", status = 1)
  }else{
    if(verbose){
        warning_message = str_glue("Parameter {name} was not specified. This script will use {name} = {default}.")
        cat(warning_message)
    }
    value = default
  }

  if(eval_string){
    value = eval(parse(text = value))
  }
  return(value)
}

# ## Help section
# In your main code, you should have something like this:
# if("--help" %in% args) {
#   cat("
 
#       Arguments:
#     --input                               - chatacter, input file
#     --out                                 - character, output file
#     --x                                   - numeric, parameter x
#     --help                                - print this text
 
#       Example:
#       Rscript this_script.R --input=input.txt --x=x --out=output.txt \n\n")
  
#   q(save="no")
# }