gglength <- function(ind, type = c("nuc", "mito"), d = "~/Projects/Botocudos/Files/Length/2017_09_15/"){
  setwd(d)
  
  #ind <- "MN00010"
  ma <- mn2ma(ind)
  #type <- "mito"
  l <- data.frame()
  for(g in type){
    tmp <- read.table(paste(ma, "_length.txt", sep = ""), header = F)
    colnames(tmp) <- c("Freq", "Length")
    tmp$Type <- g
    l <- rbind(l, tmp)
  }
  l$Type <- factor(l$Type, levels = c("mito", "nuc"))  
  colors <- c("yellowgreen", "royalblue1")
  names(colors) <- levels(l$Type)
  p <- ggplot(l, aes(x = Length, y = Freq/sum(Freq), 
                     color = Type)) +
    geom_line(size = 1.5) +
    scale_color_manual(values = colors, 
                       #breaks = c("mito", "nuc"),
                       labels = c("Mitochondrial", "Nuclear"))+
    theme(legend.position = "none") +
    labs(y = "Frequency", title = paste(ind, " (", round(boto$hits_unique_frac.endogenous.[boto$Library == ind]*100, 2),"%) ",type, sep = ""))
  return(p)
  
}
