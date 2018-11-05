library(plyr)
args = commandArgs(trailingOnly=TRUE)
panel <- read.table(args[1])
name_div <- read.table(args[2])
out <- args[3]

colnames(panel) <- c("famID", "indID", "region", "population", "label")
panel$region <- NULL
colnames(name_div) <- c("population", "region")

new_panel <- join(panel, name_div, by = "population")

print(paste("Writing new panel to", out))
print(paste("Panel contains", dim(new_panel)[1], "individuals and",
dim(new_panel)[2], "columns"))
print("Head:")
print(head(new_panel))

write.table(new_panel, out, col.names = T, row.names = F,
quote = F, sep = "\t")
