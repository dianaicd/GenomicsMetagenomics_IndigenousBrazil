#####3
# Convert tables to latex format
library(xtable)
setwd("~/Projects/Botocudos/Scripts/examples/Thesis Templates/RMarkdownTemp/sections/")
t <- read.csv("~/Desktop/SampledBotocudos_Oct312013_fromSilvia_inCopenhagenMay2017.csv")

t2 <- t[1:23,c(2,14,15,20,23,25,26)]
t2$Ethnicity <- gsub("Botocudo .*", "Botocudo", t2$Ethnicity)
colnames(t2) <- gsub(".", " ", colnames(t2))
xtable(t[1:23,c(2,13,14,15,20,23,25,26)])
print(xtable(t[1:23,c(2,13,14,15,20,23,25,26)], type = "latex"), 
      file = "sampledbotocudos.latex")
