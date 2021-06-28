args = commandArgs(trailingOnly=TRUE)
abba_path <- args[1]
corrected_path <- args[2]
output_path <- args[3]

abbababa <- read.table(abba_path, header = T)

trees <- read.table(corrected_path, header = T)[,c("H1","H2","H3")]
abbababa <- cbind(abbababa, trees)
abbababa$D <- abbababa$Numer/abbababa$Denom
abbababa <- abbababa[,c("D", "H1", "H2", "H3")]

write.table(abbababa, output_path,
col.names = T, row.names = F,
sep = "\t", quote = F)
