#library(ggplot2)
#library(cowplot)
library(scales)
ggendogenous <- function(boto, fill = "midnightblue", color = "midnightblue", lanes = F){
  rounds <- length(color)
  if(lanes){
    p <- ggplot(boto[boto$Library != "*", ],
                aes(x = Library, y = hits_unique_frac.endogenous.,
                    fill = boto[,fill])) +
      geom_hline(yintercept = seq(0.05, 0.3, 0.05), lty = "dotted",
                 size = 0.1, colour = "black") +
      geom_hline(yintercept = seq(0.075, 0.3, 0.05), lty = "dotted",
                 size = 0.1, colour = "gray") +
      geom_bar(stat = "identity", position = "dodge") +
      scale_fill_manual(values = color, name = "Sequencing") 
  }else{
    p <- ggplot(boto[boto$Library != "*", ],
                aes(x = Library, y = hits_unique_frac.endogenous.)) +
      geom_hline(yintercept = seq(0.05, 0.3, 0.05), lty = "dotted",
                 size = 0.1, colour = "black") +
      geom_hline(yintercept = seq(0.075, 0.3, 0.05), lty = "dotted",
                 size = 0.1, colour = "gray") +
      geom_bar(stat = "identity", fill = fill)
  }
  #boto <- read.table("Botocudos_summary_2017_09_18.table", header = T)
  p <- p +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_y_continuous(labels = percent(seq(0, 0.3, 0.05)), 
                       breaks = seq(0, 0.3, 0.05)) +
    labs(y = "Reads mapped to hg19\n(% of retained reads)",
         x = NULL, title = "Endogenous content") +
    geom_vline(xintercept = floor(dim(order_endo(0, 0.01, boto))[1]/rounds) + 0.5, 
               lty = "dashed", col = "gray")+
    geom_vline(xintercept = floor(dim(order_endo(0.0, 0.05, boto))[1]/rounds) + 0.5, 
               lty = "dashed", col = "gray")+
    geom_vline(xintercept = floor(dim(order_endo(0, 0.12, boto))[1]/rounds) + 0.5, 
               lty = "dashed", col = "gray")+
    geom_vline(xintercept = floor(dim(order_endo(0, 0.2, boto))[1]/rounds) + 0.5, 
               lty = "dashed", col = "gray") +
    annotate("text", x = 4, y = 0.4, label = "Endogenous content:") +
    annotate("text", x = 4, y = 0.37, label = "<1%") +
    annotate("text", x = 12, y = 0.37, label = "1% - 5%") +
    annotate("text", x = 17, y = 0.37, label = "5% - 10%") +
    annotate("text", x = 22.5, y = 0.37, label = "16% - 18%") +
    annotate("text", x = 24, y = 0.37, label = "35%")
  return(p)
  
  # hist(boto$hits_unique_frac.endogenous.[boto$Library != "*"], 
  #      breaks = c(0, 0.01, seq(0.01, 0.4, 0.02)), freq = T, xlab = "Endogenous content (fraction)", 
  #      main = "Histogram of frequencies", col = "firebrick1", border = "white") 
  # grid()

}
