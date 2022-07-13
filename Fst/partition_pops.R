library(partitions)
inds <- read.table("9ind.txt")$V1
total <- length(inds)

k = 1
for(i in 1:floor(total/2)){
  my_partitions <- setparts(c(i, total-i))
  print(dim(my_partitions))
  for(j in 1:ncol(my_partitions)){
    index_1 <- which(my_partitions[,j] == 1)
    index_2 <- which(my_partitions[,j] == 2)
    write.table(
        inds[index_1],
        paste0("partitions/group_", k, "_pop1.txt"),
        quote = F, row.names = F, col.names = F
    )
    write.table(
        inds[index_2],
        paste0("partitions/group_", k, "_pop2.txt"),
        quote = F, row.names = F, col.names = F
    )
    k = k + 1
  }
}
