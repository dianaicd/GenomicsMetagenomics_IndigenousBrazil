gglength <- function(ind, type = c("nuc", "mito"),
                     d = "~/Projects/Botocudos/Files/Length/2017_09_15/", 
                     lib = "L1",
                     gg = TRUE){
  setwd(d)
  
  #ind <- "MN00010"
  #ma <- mn2ma(ind)
  #type <- "mito"
  l <- data.frame()
  for(g in type){
    if(g == "mito"){
      tmp <- read.table(paste(d, "/",lib, ".", g, ".length", sep = ""), header = F)
    }else{
      tmp <- read.table(paste(d, "/",lib, ".length", sep = ""), header = F)
    }
    colnames(tmp) <- c("Length", "Freq")
    tmp$Type <- g
    l <- rbind(l, tmp)
  }
  l$Type <- factor(l$Type, levels = c("mito", "nuc"))  
  colors <- c("yellowgreen", "royalblue1")
  names(colors) <- levels(l$Type)
  title <- paste(ind, " (endogenous: ", 
                 round(boto$hits_unique_frac_endogenous[boto$sample == ind & 
                                                          boto$library == "All"]*100, 2),"%) ",
                 " ", type, " ", lib, sep = "")
  if(gg){
    p <- ggplot(l, aes(x = Length, y = Freq/sum(Freq), 
                       color = Type)) +
      geom_line(size = 1.5) +
      scale_color_manual(values = colors, 
                         #breaks = c("mito", "nuc"),
                         labels = c("Mitochondrial", "Nuclear"))+
      theme(legend.position = "none") +
      labs(y = "Frequency", 
           title = title)
    return(p)
  }else{
    l <- l[order(l$Length),]
    l$color <- ifelse(l$Type == "nuc", "royalblue1", "yellowgreen")
    plot(x = l$Length, y = l$Freq/sum(l$Freq), type = "l", col = l$color, bty = "n", 
         xlab = "Length", ylab = "Frequency", lwd = 1.5, 
         main = title)
    
  }
}
