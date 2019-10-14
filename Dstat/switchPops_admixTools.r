# Fix ind file for AdmixTools
library(plyr)
args = commandArgs(trailingOnly=TRUE)
panelName <- args[1]
indName <- args[2]

#indName <- "Maanasa_mask1_flip.B.P.ASM.VMAnc.ind"
#panel <- read.table("~/archive/Panels/Magic.panel", header = T)
panel <- read.table(panelName, header = T)
ind <- read.table(indName)
colnames(ind) <- c("indID", "U", "control")
ind$control <- NULL
ind$indID <- sub(":.*", "", ind$indID)
new_ind_names <- join(ind, panel, by = "indID")
new_ind_names <- new_ind_names[,c("indID", "U","population")]
dim(ind)
dim(new_ind_names)
dim(new_ind_names[complete.cases(new_ind_names),])
new_ind_names$population <- as.character(new_ind_names$population)
panel$population <- as.character(panel$population)

# this part is not great;
# if an individual is not in my Magic panel,
# it will not have a population assigned, and it will try to
# assign the population of a sample whose ID is a superstring 
# of the current individual. This is helpful for this panel
# with LS5 and LS5.variant.

# Please update your Magic panel
na <- is.na(new_ind_names$population)
myIndex <- unlist(sapply( new_ind_names$indID[na], function(x) grep(x, panel$indID)))
new_ind_names$population[na] <- panel[myIndex,]$population
dim(new_ind_names[complete.cases(new_ind_names),])
dim(ind)
new_ind_names$indID <- paste(new_ind_names$ind, new_ind_names$population, sep = ":")


write.table(new_ind_names, file=sub(".ind", ".modified.ind",indName), sep=" ", quote=F, row.names=F, col.names=F)
