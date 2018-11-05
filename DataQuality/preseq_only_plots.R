#load("~/Projects/Botocudos/Files/PreSeq/output_preseq_Sep16_all.Rda")
#boto <- read.table("~/Projects/Botocudos/Files/Summaries/Botocudos_summary_2017_09_18.table", header = T)

plot_gc_preseq <- function(ind, ylim = 1.5){
  if(ind == "MN01701"){
    blankPlot <- ggplot()+geom_blank(aes(1,1))+
      theme(
        plot.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank()
      )
    return(blankPlot)}
  
  m <- mn2ma(ind)
  f <- read.table(paste("~/Projects/Botocudos/Files/PreSeq/GC/", 
                        m,".hg19_nuc.realigned.bam.mr_gc_extrap.txt", sep = ""),
                  header = T)
  endo <- boto$hits_unique_frac.endogenous.[boto$Library == ind]
  l <- boto$hits_length.nuclear.[boto$Library == ind]
  sequenced <- boto$seq_reads_se[boto$Library==ind]
  txt <- paste("Sequenced:", prettyNum(sequenced, big.mark = ","), sep = "\n")
  
  p <- ggplot(f,
              aes(x = (TOTAL_BASES/l)/
                    (boto$hits_unique.nuclear.[
                      boto$Library==ind]/
                       boto$seq_reads_se[boto$Library==ind]), 
                  y = (EXPECTED_COVERED_BASES)/(1e+9)/2.7,
                  ymin = (EXPECTED_COVERED_BASES - LOWER_95.CI)/(1e+9)/2.7,
                  ymax = (EXPECTED_COVERED_BASES + LOWER_95.CI)/(1e+9)/2.7)) +
    geom_line() +
    geom_ribbon(alpha = 0.6, fill = "salmon")  +
    coord_cartesian(xlim = c(0, 3e+8), ylim = c(0, ylim)) +
    labs(title = paste(ind, " (", 
                       percent(signif(boto$hits_unique_frac.endogenous.[
                         boto$Library == ind], 2)), ")", sep = ""),
         x = "Lanes", y = "Fraction of genome") +
    theme(legend.position = "none") +
    annotate("text", x = 9e+7, y = 0.9*ylim, label = txt)
  
  return(p)
}

plot_expreads_preseq <- function(ind, ylim, verify = F, xlim = 3e8, useboto = T){
  
  # if(ind == "MN01701"){
  #   blankPlot <- ggplot()+geom_blank(aes(1,1))+
  #     theme(
  #       plot.background = element_blank(), 
  #       panel.grid.major = element_blank(),
  #       panel.grid.minor = element_blank(), 
  #       panel.border = element_blank(),
  #       panel.background = element_blank(),
  #       axis.title.x = element_blank(),
  #       axis.title.y = element_blank(),
  #       axis.text.x = element_blank(), 
  #       axis.text.y = element_blank(),
  #       axis.ticks = element_blank(),
  #       axis.line = element_blank()
  #     )
  #   return(blankPlot)}
  
  # Be careful and remember whether you are using a MA ID or a MN ID
  if(!is.data.frame(useboto)){
    sequenced <- boto$seq_reads_se[boto$Target==ind]
    conv_factor <- sequenced/
      boto$seq_retained_reads[boto$Target == ind]
    useboto <- boto
  }else{
    sequenced <- useboto$seq_reads_se[useboto$Target==ind]
    conv_factor <- sequenced/
      useboto$seq_retained_reads[useboto$Target == ind]
  }


  
  txt <- paste("Sequenced:", prettyNum(sequenced, big.mark = ","), sep = "\n")
  p <- ggplot(semi_final[[ind]]$Reads, 
              aes(x = Total_reads*conv_factor, 
                  y = Expected,
                  ymin = Expected - sqrt(Var),
                  ymax = Expected + sqrt(Var))) +
    geom_line() +
    coord_cartesian(xlim = c(0, xlim), ylim = c(0, ylim)) +
    geom_ribbon(alpha = 0.6, fill = "chartreuse3")  +
    geom_vline(xintercept = useboto$seq_reads_se[useboto$Library==ind], linetype = "dashed", color = "gray") +
    labs(title = paste(ma2mn(ind), " (", 
                       percent(signif(boto$hits_unique_frac.endogenous.[
                         boto$Target == ind], 2)), ")", sep = ""), 
         x = "Reads", y = "Expected unique reads") +
    theme(legend.position = "none") +
    annotate("text", x = 0.15*xlim, y = 0.9*ylim, label = txt)
  
  if(verify){
    #histo <- read.table(paste("~/Projects/Botocudos/Files/Histograms/2018_Sam/", mn2ma(ind), ".bam_duphist.txt", sep = ""))
    histo <- read.table(paste("~/Projects/Botocudos/Files/Histograms/2018_07_23/", ind, "_histogram.txt", sep = ""))
    p <- p + geom_point(x = useboto$seq_reads_se[useboto$Target==ind],
                        y = sum(histo$V2), color = "red", size = 2)
  }
  return(p)
}
