# Red: C to T substitutions.
# Blue: G to A substitutions.
# Grey: All other substitutions.
# Orange: Soft-clipped bases.
# Green: Deletions relative to the reference.
# Purple: Insertions relative to the reference.


##### Plot
damage_plot <- function(ind, type, lib, d ="~/Projects/Botocudos/Files/mapDamage/2018_07_23/",
                        asis = F, sufix = "_Human_results_mapDamage", plot.title  = F){
  mean_muts <- function(mutation, ref){
    mutation <- as.character(mutation)
    ref <- as.character(ref)
    mean_pos <- function(pos){
      res <- (md[md$End == pos & md$Std == "+", mutation]/
                ifelse(md[md$End == pos & md$Std == "+", ref] > 0, md[md$End == pos & md$Std == "+", ref], 1 ) +
                md[md$End == pos & md$Std == "-", mutation]/
                ifelse(md[md$End == pos & md$Std == "-", ref] >0, md[md$End == pos & md$Std == "-", ref], 1 ))/2
      return(res)
    }
    mut_5p <- mean_pos("5p")
    mut_3p <- mean_pos("3p")
    return(data.frame("left" = mut_5p, "right" = mut_3p))
  }
  
  #setwd("~/Projects/Botocudos/Files/mapDamage/2017_09_15/")
  setwd(d)
  if(asis){
    ma <- ind
  }else{
    ma <- mn2ma(ind)
  }

  # md <- read.table(paste(ma, ".hg19_", type, 
  #                        ".mapDamage/", ind, "/misincorporation.txt", sep = ""),
  #                  header = T)
  # 
   md <- read.table(paste(ma, sufix,
                          "/misincorporation.txt", sep = ""),
                    header = T)
  mean_GA <-  mean_muts("G.A", "G")
  #12 - 21
  
  mean_CT <-  mean_muts("C.T", "C")
  
  # Soft clipped
  mean_soft <- mean_muts("S", "Total")
  
  # All other mutations
  other <- data.frame(mut = colnames(md)[12:21], ref = gsub("\\..*", "", colnames(md)[12:21]))
  tmp <- lapply(seq(1, 10), function(i) mean_muts(other$mut[i], other$ref[i]))
  
  mean_other <- tmp[[1]]
  
  for(i in seq(1,10)){
    
    mean_other$right <- rowMeans(cbind(mean_other$right, tmp[[i]]$right), na.rm = T)
    mean_other$left <- rowMeans(cbind(mean_other$left, tmp$left[[i]]), na.rm = T)
  }
  
  # Deletions
  del <- data.frame(mut = colnames(md)[22:25], ref = rep("Total",4))
  tmp <- lapply(seq(1, 10), function(i) mean_muts(del$mut[i], del$ref[i]))
  
  mean_del <- tmp[[1]]
  
  for(i in seq(1,10)){
    
    mean_del$right <- rowMeans(cbind(mean_del$right, tmp[[i]]$right), na.rm = T)
    mean_del$left <- rowMeans(cbind(mean_del$left, tmp$left[[i]]), na.rm = T)
  }
  
  # Insertions
  ins <- data.frame(mut = colnames(md)[26:29], ref = rep("Total", 4))
  tmp <- lapply(seq(1, 10), function(i) mean_muts(ins$mut[i], ins$ref[i]))
  
  mean_ins <- tmp[[1]]
  
  for(i in seq(1,10)){
    
    mean_ins$right <- rowMeans(cbind(mean_ins$right, tmp[[i]]$right), na.rm = T)
    mean_ins$left <- rowMeans(cbind(mean_ins$left, tmp$left[[i]]), na.rm = T)
  }
  if(plot.title == F){
    plot.title <-  paste(ind, lib, sep = ", ")
  }
  #####
#  par(mfrow=c(1,2))
  par(mar = c(4, 3, 3, 1))
  plot(x = seq(1,25), 
       y = mean_GA$left[1:25], type = "l", xlim = c(1,25), ylim = c(0, 0.3), 
       col = NA, lwd = 3, axes = F, ylab = "Frequency", xlab = "", 
       main = plot.title)
  axis(side = 1, at = seq(0,25), las = 2)
  axis(side = 2, las =2)
  lines(x = seq(1,25), 
        y = mean_other$left[1:25], xlim = c(0, 25), ylim = c(0, 0.3), 
        col = "grey", bty = "n", lwd = 1)
  lines(x = seq(1,25), 
        y = mean_soft$left[1:25], xlim = c(0, 25), ylim = c(0, 0.3), 
        col = "orange", bty = "n", lwd = 1)
  lines(x = seq(1,25), 
        y = mean_del$left[1:25], xlim = c(0, 25), ylim = c(0, 0.3), 
        col = "green", bty = "n", lwd = 1)
  lines(x = seq(1,25), 
        y = mean_ins$left[1:25], xlim = c(0, 25), ylim = c(0, 0.3), 
        col = "purple", bty = "n", lwd = 1)
  lines(x = seq(1,25), 
        y = mean_CT$left[1:25], xlim = c(0, 25), ylim = c(0, 0.3), 
        col = "red", bty = "n", lwd = 3)
  lines(x = seq(1,25), 
        y = mean_GA$left[1:25], type = "l", xlim = c(1,25), ylim = c(0, 0.3), 
        col = "blue", lwd = 3)
  
  par(mar = c(4,1,3,3))
  plot(x = seq(-25, -1), 
       y = mean_GA$right[25:1], type = "l", xlim = c(-25, 0), ylim = c(0, 0.3), 
       col = "blue", bty = "c", lwd = 3, ylab = "", xlab = "",
       axes = F, main = "")
  axis(side = 1, at = seq(-25, -1), las = 2)
  axis(side = 4, las =2)
  lines(x = seq(-25,-1), 
        y = mean_other$left[1:25], xlim = c(-25,0), ylim = c(0, 0.3), 
        col = "grey", bty = "n", lwd = 1)
  lines(x = seq(-25,-1), 
        y = mean_soft$left[1:25], xlim = c(-25,0), ylim = c(0, 0.3), 
        col = "orange", bty = "n", lwd = 1)
  lines(x = seq(-25,-1), 
        y = mean_del$left[1:25], xlim = c(-25,0), ylim = c(0, 0.3), 
        col = "green", bty = "n", lwd = 1)
  lines(x = seq(-25,-1), 
        y = mean_ins$left[1:25], xlim = c(-25,0), ylim = c(0, 0.3), 
        col = "purple", bty = "n", lwd = 1)
  lines(x = seq(-25, -1), 
        y = mean_CT$right[25:1], type = "l", xlim = c(-25, 0), ylim = c(0, 0.3), 
        col = "red", bty = "n", lwd = 3)
  
}

#==============================================================================#
# Plot with input from bamdamage

plot_bamdamage <- function(five_prime, three_prime, cex = 1, 
                           y2 = 0.35, plot.title, myColor = c(),
                           legend_only = F, ncol_legend = 3, legend_cex = 1){
  
  if(length(myColor) == 0){
    myColor <- c("royalblue", "firebrick3", "lightblue3", "indianred1",
                 "orange","orange3",
                 "plum", "plum4",
                 "paleturquoise3", "paleturquoise4",
                 "rosybrown2", "rosybrown3")
    names(myColor) <- c("C..T", "G..A","T..C", "A..G",
                        "C..A", "C..G",
                        "T..A", "T..G", 
                        "A..C", "A..T",
                        "G..C", "G..T")
  }
  
  if(legend_only){
    myIndex <- round(seq(1, 13, length.out = ncol_legend+1))
    for(i in 1:(length(myIndex)-1)){
      plot(1:10,1:10, type = "n", bty = "n", axes = F, xlab = NA, ylab = NA)
      
      if(i == 1){
        legend(1,10, legend = sub("\\.\\.", " to ",
                                  names(myColor)[myIndex[i]:(myIndex[i+1]-1)]),
               col = myColor[myIndex[i]:(myIndex[i+1]-1)],
               lty = 1, bty = "n", cex = legend_cex, lwd = 1.6,
               title = "Substitution type")
      }else{
        legend(1,10, legend = sub("\\.\\.", " to ",names(myColor)[myIndex[i]:(myIndex[i+1]-1)]),
               col = myColor[myIndex[i]:(myIndex[i+1]-1)],
               lty = 1, bty = "n", cex = legend_cex, lwd = 1.6
        )
      }
    }
    return()
  }
  
  plot_lines <- function(damage, left = T){
    plot(damage$C..T[1:25], type = "n",
         bty = "n", 
        ylab = "Frequency", xlab = NA,
         cex = cex, axes = F, ylim = c(0, y2))
    if(left){
      title(main = plot.title, xlab = "5'-position")
      
      axis(side = 1, at = seq(0,25,5))
      axis(side=2)
    }else{
      axis(side = 1, at = seq(1,25,5), labels = seq(25, 1,-5))
      axis(side = 4)
      title(xlab="3'-position")
    }
    
    for(type in names(myColor)){
      lines(damage[1:25, type],
            col = alpha(myColor[type], 0.8), 
            lwd = 1.5, 
            cex = cex)
    }
  }
  
  par(mar = c(4, 4, 3, 1))
  plot_lines(damage = five_prime)
  
  par(mar = c(4,1,3,3))
  plot_lines(damage = three_prime[25:1,], left = F)
}