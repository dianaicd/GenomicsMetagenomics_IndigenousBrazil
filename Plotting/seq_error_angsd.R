# Plot error rates from angsd output

#------------------------------------------------------------------------------#
# a trick I'll use a few times
repeat_aes <- function(myVec, nInd, mix = F)  {
  if(mix){
    myNewVec <- rep(myVec, round(nInd/length(myVec)), length.out = nInd)
  }else{
    myNewVec <- rep(myVec, each = round(nInd/length(myVec)), length.out = nInd)
  }
  return(myNewVec)
}

#------------------------------------------------------------------------------#
# function to return a data frame with the specification for the legend
legend_error <- function(sm, myColors = c(), nCol = 6, add_texture = F, 
                         myDensities =  c(10,20,30,7),
                         myAngles =  c(0,45,-45,36)){
  mySamples <- read.table(sm)
  nInd <- nrow(mySamples)
  if(length(myColors) == 0){
    require("RColorBrewer")
    myColors <- brewer.pal(nCol, "Dark2")  
  }
  all_colors <- repeat_aes(myColors, nInd)
  if(add_texture){
    all_densities <- repeat_aes(myDensities, nInd = nInd, mix = T)
    all_angles <- repeat_aes(myAngles, mix = T, nInd = nInd)
    myLegend <- data.frame(sample = mySamples$V1,
                           color = all_colors,
                           density = all_densities,
                           angle = all_angles)
  }else{
    myLegend <- data.frame(sample = mySamples$V1,
                           color = all_colors)
  }
  return(myLegend)
}

#------------------------------------------------------------------------------#
# function to plot the type-specific error rates
plot_bars_error <- function(error, sm, myColors = c(), nCol = 6,
                            add_texture = F,
                            myDensities =  c(10,20,30,7),
                            myAngles =  c(0,45,-45,36), 
                            myMain = "", return_legend = F,
                            ymax = F){
  mySamples <- read.table(sm)
  nInd <- nrow(mySamples)
  dat <- read.table(error, skip = 0, header = TRUE, nrows = nInd)

  if(nCol*length(unique(myAngles)) < nInd){
    message("-------------- Not enough colors and/or angles!!! --------------")
    break()
  }
  
  if(!ymax){
    ymax <- max(dat)  
  }
  
  names <- colnames(dat)

  labels = c("C->A | G->T", "G->A | C->T", 
             "T->A | A->T\nA->C | T->G\nG->C | C->G\nT->C | A->G")
  
  if(length(myColors) == 0){
    require("RColorBrewer")
    myColors <- brewer.pal(nCol, "Dark2")  
  }
  
  all_colors <- repeat_aes(myColors, nInd)

  # Plot solid bars
  barplot(as.matrix(dat[,1:3]), beside=T,
           ylim=c(0,ymax), col = all_colors,
          cex.names=0.6, main = myMain, 
          ylab="Error rate",  names.arg=rep("",3),
          border = NA, add = F)
  
  # Add texture (white lines)
  if(add_texture){
    all_densities <- repeat_aes(myDensities, nInd = nInd, mix = T)
    all_angles <- repeat_aes(myAngles, mix = T, nInd = nInd)
    
    barplot(as.matrix(dat[,1:3]),beside=T,
            ylim=c(0,ymax),
            col = "white", density = all_densities,
            angle = all_angles, names.arg=rep("",3), add = T, 
            border = "white")
    
    myLegend <- data.frame(sample = mySamples$V1,
                           color = all_colors,
                           density = all_densities,
                           angle = all_angles)
  }else{
    myLegend <- data.frame(sample = mySamples$V1,
                           all_colors)
  }
  
  grid(nx=0, ny=NULL)
  
  pos <- seq(round(nInd/2), by = nInd, length.out = 3)
  mtext(labels, 1, at = pos, cex=0.8, padj = 1)
  
}

#------------------------------------------------------------------------------#
# function to plot average error rates

plot_avg_error <- function(error, nInd, myColors = "gray"){
  dat <- read.table(error, skip = nInd+1, header = F, sep = "\t")
  
  myPercent <- data.frame(t(unlist(apply(dat, 1, function(x) unlist(strsplit(x, " "))))))
  myPercent[,2] <- as.numeric(as.character(myPercent[,2]))
  rownames(myPercent) <- myPercent$V11
  
  barplot(myPercent[,2], beside = T, las = 2, width = 1, 
          main = "Average error rate", col = myColors)
  
  sp <- 1.2
  text(x = seq(sp/2, nInd*sp, sp), par("usr")[3] - 0.01,
       labels = myPercent[,1], srt = 90, pos = 1, xpd = TRUE)
}




