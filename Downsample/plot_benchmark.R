library(ggplot2)

pipelines <- c("me", "fred", "shuf")
reads <- c("100","1000", "10000", "100000", "1000000", "10000000", "100000000")

path <- "~/Projects/Botocudos/Files/Downsample/2019_08_13/benchmarks/"

myRes <- data.frame()
for(n in reads){

  for(p in pipelines){
    new_path <- paste(path, p, "_MN0008_", n, ".benchmark.txt", sep = "")
    tmp_df <- read.table(new_path, header = T)
    tmp_df$pipeline <- p
    tmp_df$reads <- n
    myRes <- rbind(myRes, tmp_df[,c("s", "pipeline", "reads")])
  }
  
  
}

ggplot(myRes, aes(reads, s, col = pipeline)) +
  geom_boxplot()
