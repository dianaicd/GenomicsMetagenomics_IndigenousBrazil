library(ggplot2)
library(cowplot)
source("~/Projects/Botocudos/Scripts/translate_ids.R")
relat <- read.table("~/Projects/Botocudos/Files/Relatedness/2017_09_15/gl.res", header = T)
bam <- read.table("~/Projects/Botocudos/Files/Relatedness/2017_09_15/bam.list")
bam$V1 <- gsub(".*/", "", bam$V1)
bam$V1 <- gsub("\\..*", "", bam$V1)
bam$V1 <- sapply(bam$V1, function(x) ma2mn(x))

ggplot(relat, aes(x = a, y = b)) +
  geom_tile(aes(fill = k0)) +
  geom_tile(aes(x = b, y = a, fill = k0)) +
  scale_fill_gradient2(low="white", high="darkred", guide="colorbar") +
  scale_x_discrete(limits = seq(0,22), labels = bam$V1) +
  theme(axis.text.x = element_text(angle = 270)) +
  labs(x = NULL, y = NULL, title = expression(k[0])) +
  scale_y_discrete(limits = seq(0,22), labels = bam$V1) 

ggplot(relat, aes(x = a, y = b)) +
  geom_tile(aes(fill = k1)) +
  geom_tile(aes(x = b, y = a, fill = k1)) +
  scale_fill_gradient2(low="white", high="darkred", guide="colorbar") +
  scale_x_discrete(limits = seq(0,22), labels = bam$V1) +
  theme(axis.text.x = element_text(angle = 270)) +
  labs(x = NULL, y = NULL, title = expression(k[1])) +
  scale_y_discrete(limits = seq(0,22), labels = bam$V1) 

ggplot(relat, aes(x = a, y = b)) +
  geom_tile(aes(fill = k2))  +
  geom_tile(aes(x = b, y = a, fill = k2)) +
  scale_fill_gradient2(low="white", high="darkred", guide="colorbar") +
  scale_x_discrete(limits = seq(0,22), labels = bam$V1) +
  theme(axis.text.x = element_text(angle = 270)) +
  labs(x = NULL, y = NULL, title = expression(k[2])) +
  scale_y_discrete(limits = seq(0,22), labels = bam$V1) 
