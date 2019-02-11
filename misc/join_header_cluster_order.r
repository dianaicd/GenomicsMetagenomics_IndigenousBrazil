library(plyr)
file <- "Maanasa_mask1_flip.americas.names"
file2 <- "Maanasa_mask1_flip.americas.clust"
m <- "Maanasa.americas.Moreno_names.txt"
h <- read.table(file)
colnames(h) <- c("id")
clus <- read.table(file2)
colnames(clus) <- c("id", "id2", "pop")
clus$id <- paste(clus$id, clus$id, sep = "_")

moreno <- read.table(m)
colnames(moreno) <- c("pop", "region")

head_clus <- join(h, clus, by = "id")
head(head_clus)

panel <- join(head_clus, moreno, by = "pop")
head(panel)

write.table(panel, "America_Maanasa_ordered.txt", sep = "\t", col.names = F, row.names = F, quote = F)