# Fix ind file for AdmixTools
library(plyr)
args = commandArgs(trailingOnly=TRUE)
panel <- args[1]
indName <- args[2]

#indName <- "Maanasa_mask1_flip.B.P.ASM.VMAnc.ind"
#panel <- read.table("~/archive/Panels/Magic.panel", header = T)
ind <- read.table(indName)
colnames(ind) <- c("indID", "U", "control")
ind$control <- NULL
ind$indID <- sub(":.*", "", ind$indID)
x <- join(ind, panel, by = "indID")
x <- x[,c("indID", "U","population")]
dim(ind)
dim(x)
dim(x[complete.cases(x),])
x$population <- as.character(x$population)
panel$population <- as.character(panel$population)
x$population[is.na(x$population)] <- panel[unlist(sapply( x$indID[is.na(x$population)], 
function(x) grep(x, panel$indID))),]$population
dim(x[complete.cases(x),])
dim(ind)
x$indID <- paste(x$ind, x$population, sep = ":")


write.table(x, file=sub(".ind", ".modified.ind",indName), sep=" ", quote=F, row.names=F, col.names=F)
