library(plyr)
args = commandArgs(trailingOnly=TRUE)

# Master table
master <- read.table(args[1], header = T)
# IDs of the individuals in this analysis
ids <- read.table(args[2], header = T)

merged <- join(ids, master, by = "ID")

complete <- dim(merged[complete.cases(merged),])
if(dim(merged)[1] != complete[1]){
    missing <- dim(merged)[1] - complete[1]
    print(paste("Warning: ", missing, "missing individuals in the master table."))
}

print(merged$ID[!complete.cases(merged)])

write.table(merged, paste("panel", args[2], sep = "_"), sep = "\t",
col.names = T, row.names = F, quote = F)