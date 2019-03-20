args = commandArgs(trailingOnly=TRUE)

preffix <- args[1]
panel_name <- args[2]
ind_name <- args[3]
n <- args[4]
program <- args[5]

suffix <- paste( panel_name, ind_name, sep = "_")
files <- list.files(".", pattern=paste(preffix, suffix, sep = ".*"))
length(files)
toscp <- c()
ourbest <- c()
for(k in seq(2, n)){
  file_name <- paste(preffix, "_k", k, "_", panel_name, "_", ind_name, sep = "")
  file <- read.csv(file_name, sep = "\t")
  like <- max(file$V2)
  index <- which(file$V2 == like)
  print(paste("Best likelihood for k", k, "is", like, "of replicate", file$V1[index]))

  if(program == "ADMIXTURE"){
    best <- paste(k, "/", panel_name, "_", sub(".txt", "",ind_name),
                  "_",file$V1[index], ".", k, sep = "")
  }else{
    best <- paste(k, "/", panel_name, "_", sub(".txt", "",ind_name),
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
  print(paste("scp dcruzdva@prd.vital-it.ch:", curdir, "/{",
        toscp, "}.Q ./", sep = ""))
}else{
  print(paste("scp dcruzdva@prd.vital-it.ch:", curdir, "/{",
        toscp, "}.qopt ./", sep = ""))
}

if(program == "ADMIXTURE"){
  forloop <- paste("for(k in paste(c(",  ourbest,  "), seq(2,",  n, "), sep = '.'))"  )
}else{
  forloop <- paste("for(k in paste(seq(2,",  n,  "), c(",  ourbest,  "), sep = '_'))"  )
}


print(forloop)
for(k in paste(seq(2,n), ourbest, sep = "_")){
  print(k)
}
