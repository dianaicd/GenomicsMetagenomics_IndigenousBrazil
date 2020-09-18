args = commandArgs(trailingOnly=TRUE)

prefix <- args[1]
panel_name <- args[2]
# ind_name <- args[3]
min_k <- args[3]
max_k <- args[4]
program <- args[5]

# sufix <- paste( panel_name, ind_name, sep = "_")
files <- list.files(".", pattern = paste(panel_name, sep = ".*"))
length(files)
toscp <- c()
ourbest <- c()
for(k in seq(min_k, max_k)){
  file_name <- paste(prefix, "_k", k, "_", panel_name, ".txt",#"_", ind_name, 
  sep = "")
  file <- read.csv(file_name, sep = "\t", header = F)
  like <- max(file$V2)
  index <- which(file$V2 == like)
  print(paste("Best likelihood for k", k, "is", like, "of replicate", file$V1[index]))

  if(program == "ADMIXTURE"){
    best <- paste(k, "/", panel_name, #"_", sub(".txt", "",ind_name),
                  "_",file$V1[index], ".", k, sep = "")
  }else{
    best <- paste(k, "/", panel_name, #"_", sub(".txt", "",ind_name),
                  "_k", k, "_",file$V1[index], sep = "")
  }

  toscp <- paste(toscp, best, sep = ",")
  ourbest <- paste(ourbest, file$V1[index], sep = ",")
}
#   for(k in paste(seq(2,10), c(201,1,252,50,231,257,248,227,56), sep = "_"))
ourbest <- sub(",", "", ourbest)

toscp <- sub(",", "", toscp)
curdir <- getwd()
if(program == "ADMIXTURE"){
  print(paste("scp Axiom:", #dcruzdva@prd.vital-it.ch:", 
  curdir, "/{",
        toscp, "}.Q ./", sep = ""))
}else{
  scp_command <- paste("scp Axiom:",  curdir, "/{", toscp, "}.qopt ./", sep = "")
  mv_commands <- sapply( sub(".*/", "", unlist( strsplit(toscp, ",") ) ), 
                        function( x ) 
                                            paste( "mv ", x, ".qopt ", 
                                                    sub("_\\d+$", ".qopt", x), 
                                            sep = "" ) )
  print(scp_command)
  for( command in mv_commands ){ print(command) }
}

if(program == "ADMIXTURE"){
  forloop <- paste("for(k in paste(c(",  ourbest,  "), seq(2,",  max_k, "), sep = '.'))"  )
}else{
  forloop <- paste("for(k in paste(seq(",min_k,",",  max_k,  "), c(",  ourbest,  "), sep = '_'))"  )
}


# print(forloop)
# for(k in paste(seq(2, max_k), ourbest, sep = "_")){
#   print(k)
# }
