plot_bases_covered <- function(boto){
  
  m <- max(boto$Bases_covered)
  ggplot(boto[boto$Library != "*",], 
         aes(x = Library, y = Bases_covered)) +
    geom_hline(yintercept = seq(5e7, 1.5e8, 5e7), 
               lty = "dotted", size = 0.1, colour = "black") +
    geom_hline(yintercept = seq(7.5e7, 1.5e8, 5e7), 
               lty = "dotted", size = 0.1, colour = "gray") +
    geom_bar(stat = "identity", fill = "lightcoral") +
    theme(axis.text.x = element_text(angle = 90)) +
    labs(y = "Number of unique bases", title = "Bases covered", x = NULL) +
    geom_vline(xintercept = dim(order_endo(0, 0.01, boto))[1] + 0.5, 
               lty = "dashed", col = "gray")+
    geom_vline(xintercept = dim(order_endo(0.0, 0.05, boto))[1] + 0.5, 
               lty = "dashed", col = "gray")+
    geom_vline(xintercept = dim(order_endo(0, 0.12, boto))[1] + 0.5, 
               lty = "dashed", col = "gray")+
    geom_vline(xintercept = dim(order_endo(0, 0.2, boto))[1] + 0.5, 
               lty = "dashed", col = "gray") +
    annotate("text", x = 4, y = 1.2*m, label = "Endogenous content:") +
    annotate("text", x = 4, y = 1.1*m, label = "<1%") +
    annotate("text", x = 10, y = 1.1*m, label = "1% - 5%") +
    annotate("text", x = 16, y = 1.1*m, label = "5% - 10%") +
    annotate("text", x = 21.5, y = 1.1*m, label = "16% - 18%") +
    annotate("text", x = 23, y = 1.1*m, label = "35%") 
}

plot_depth_cov <- function(boto){
  m <- max(boto$hits_coverage.nuclear.)
  ggplot(boto[boto$Library != "*",], 
         aes(x = Library, y = hits_coverage.nuclear.)) +
    geom_hline(yintercept = seq(0.02, 0.06, 0.02), lty = "dotted", colour = "black", size = 0.1) +
    geom_hline(yintercept = seq(0.03, 0.06, 0.02), lty = "dotted", colour = "gray", size = 0.1) +
    geom_bar(stat = "identity", aes(fill = "hits_coverage.nuclear.")) +
    geom_bar(aes(y = Bases_covered_fraction, fill = "Bases_covered_fraction"), 
             stat = "identity") +
    labs(y = NULL, title = "Depth of coverage", x = NULL) +
    scale_fill_manual(name = "", 
                      breaks = c("hits_coverage.nuclear.",
                                 "Bases_covered_fraction"),
                      values = c("hits_coverage.nuclear."="black", 
                                 "Bases_covered_fraction"="mediumvioletred"),
                      labels = c("Depth of coverage", 
                                 "Fraction of the genome covered")) +
    geom_vline(xintercept = dim(order_endo(0, 0.01, boto))[1] + 0.5, 
               lty = "longdash", col = "gray")+
    geom_vline(xintercept = dim(order_endo(0.0, 0.05, boto))[1] + 0.5, 
               lty = "longdash", col = "gray")+
    geom_vline(xintercept = dim(order_endo(0, 0.12, boto))[1] + 0.5, 
               lty = "longdash", col = "gray")+
    geom_vline(xintercept = dim(order_endo(0, 0.2, boto))[1] + 0.5, 
               lty = "longdash", col = "gray") +
    annotate("text", x = 4, y = 1.2*m, label = "Endogenous content:") +
    annotate("text", x = 4, y = 1.1*m, label = "<1%") +
    annotate("text", x = 10, y = 1.1*m, label = "1% - 5%") +
    annotate("text", x = 16, y = 1.1*m, label = "5% - 10%") +
    annotate("text", x = 21.5, y = 1.1*m, label = "16% - 18%") +
    annotate("text", x = 23, y = 1.1*m, label = "35%")  +
    theme(axis.text.x = element_text(angle = 90), legend.position = "top")
}
