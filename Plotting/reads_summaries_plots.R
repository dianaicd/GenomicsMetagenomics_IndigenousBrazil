#library(ggplot2)
#library(cowplot)

summary_reads_plot <- function(boto){
  m <- max(boto$seq_reads_se)

  ggplot(boto[boto$Library != "*", ], aes(x = Library)) +
    geom_bar(aes(y = seq_reads_se,
                 fill = "seq_reads_se"), stat = "identity") +
    geom_bar(aes(y = seq_retained_reads, fill = "seq_retained_reads"), 
             stat = "identity") +
    geom_bar(aes(y = hits_unique_nuclear, fill = "hits_unique_nuclear"), 
             stat = "identity") +
    labs(x=NULL, y="Number of reads", title = "Summary of sequenced reads") +
    scale_fill_manual(name = "Reads:", 
                      breaks = c("seq_reads_se", "seq_retained_reads", 
                                 "hits_unique_nuclear"),
                      values = c("seq_reads_se"="gray", 
                                 "seq_retained_reads"="cadetblue", 
                                 "hits_unique_nuclear"="chocolate"),
                      labels = c("Sequenced", "Retained", "Aligned to Human genome")) +
    theme(axis.text.x = element_text(angle = 90), legend.position = "top") +
    geom_vline(xintercept = dim(order_endo(0, 0.01, boto))[1] + 0.5, 
               lty = "dashed", col = "gray")+
    geom_vline(xintercept = dim(order_endo(0.0, 0.05, boto))[1] + 0.5, 
               lty = "dashed", col = "gray")+
    geom_vline(xintercept = dim(order_endo(0, 0.12, boto))[1] + 0.5, 
               lty = "dashed", col = "gray")+
    geom_vline(xintercept = dim(order_endo(0, 0.2, boto))[1] + 0.5, 
               lty = "dashed", col = "gray") +
    geom_vline(xintercept = dim(order_endo(0, 0.4, boto))[1] + 0.5, 
               lty = "dashed", col = "gray") +
    geom_vline(xintercept = dim(order_endo(0, 0.6, boto))[1] + 0.5, 
               lty = "dashed", col = "gray") +
    annotate("text", x = 4, y = 1.2*m,#4e+7, 
             label = "Endogenous content:") +
    annotate("text", x = 4, y = 1.1*m, #3.5e+7, 
             label = "<1%") +
    annotate("text", x = 12, y = 1.1*m,#3.5e+7,
             label = "1% - 5%") +
    annotate("text", x = 17, y = 1.1*m,#3.5e+7,
             label = "5% - 10%") +
    annotate("text", x = 21.5, y = 1.1*m,#3.5e+7,
             label = "16% - 18%") +
    annotate("text", x = 23, y = 1.1*m,#3.5e+7,
             label = "35%") +
    annotate("text", x = 24, y = 1.1*m,#3.5e+7,
             label = "50%")

}

hg_mapping_plot <-function(boto, ylim = F){
  if(!is.numeric(ylim)){
    m <- max(boto$hits_unique_nuclear)
    ylim <- 1.3*m
  }else{
    m <- ylim/1.3
  }

  ggplot(boto[boto$Library != "*", ], aes(x = Library,
                                        y = hits_unique_nuclear)) +
  geom_bar(stat = "identity", fill = "chocolate") +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x=NULL, y="Mapped reads", title = "Human genome") +
    geom_vline(xintercept = dim(order_endo(0, 0.01, boto))[1] + 0.5, 
               lty = "dashed", col = "gray")+
    geom_vline(xintercept = dim(order_endo(0.0, 0.05, boto))[1] + 0.5, 
               lty = "dashed", col = "gray")+
    geom_vline(xintercept = dim(order_endo(0, 0.12, boto))[1] + 0.5, 
               lty = "dashed", col = "gray")+
    geom_vline(xintercept = dim(order_endo(0, 0.2, boto))[1] + 0.5, 
               lty = "dashed", col = "gray")+
    annotate("text", x = 4, y = 1.2*m, label = "Endogenous content:") +
    annotate("text", x = 4, y = 1.1*m, label = "<1%") +
    annotate("text", x = 10, y = 1.1*m, label = "1% - 5%") +
    annotate("text", x = 16, y = 1.1*m, label = "5% - 10%") +
    annotate("text", x = 21.5, y = 1.1*m, label = "16% - 18%") +
    annotate("text", x = 23, y = 1.1*m, label = "35%") +
    coord_cartesian(ylim = c(0, ylim))
}
# ggplot(boto[boto$Library != "*", ], aes(x = Library,
#                                         y = hits_unique_frac_nuclear)) +
#   geom_bar(stat = "identity") +
#   theme(axis.text.x = element_text(angle = 90)) +
#   labs(x=NULL, y="Reads mapping to Human genome\n(fraction of retained)")
