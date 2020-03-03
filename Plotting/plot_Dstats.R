# Functions to make plots related to D-stats
switch_h1_h2 <- function(Dstat, pop_in_h1){
  idx_to_switch <- which(Dstat$H2 == pop_in_h1)
  if(length(idx_to_switch)){
    h2 <- Dstat$H1[idx_to_switch]
    Dstat$H1[idx_to_switch] <- pop_in_h1
    Dstat$H2[idx_to_switch] <- h2
    Dstat$Z[idx_to_switch] <- -Dstat$Z[idx_to_switch]
    Dstat$D[idx_to_switch] <- -Dstat$D[idx_to_switch]
    if("JK.D" %in% colnames(Dstat)){
      Dstat$JK.D[idx_to_switch] <- -Dstat$JK.D[idx_to_switch]  
    }
    
  }
  return(Dstat)
}

switch_h2_h1 <- function(Dstat, pop_in_h2){
  idx_to_switch <- which(Dstat$H1 == pop_in_h2)
  if(length(idx_to_switch)){
    h1 <- Dstat$H2[idx_to_switch]
    Dstat$H2[idx_to_switch] <- pop_in_h2
    Dstat$H1[idx_to_switch] <- h1
    
    Dstat$Z[idx_to_switch] <- -Dstat$Z[idx_to_switch]
    Dstat$D[idx_to_switch] <- -Dstat$D[idx_to_switch]
    if("JK.D" %in% colnames(Dstat)){
      Dstat$JK.D[idx_to_switch] <- -Dstat$JK.D[idx_to_switch]  
    }
  }
  return(Dstat)
}

moreno_S35A <- function(abba, pops_in_title = c("SpCave", "Mixe", "Surui", "LagoaSta"), 
                        pops_in_y = c("USR1", "Anzick1", "Mixe", "Karitiana",
                                                 "Surui", "LagoaSanta", "MN0008", 
                                                 "Aconcagua", "AncKaweskar", "AncYamana", 
                                                 "Ayayema", "Aymara", "Chane", "TrailCreek", 
                                                 "Huichol", "Lovelock", "SouthWestOntario", 
                                                 "Maya", "Piapoco", "Pima", "PuntaSantaAna",
                                                 "Quechua", "SpiritCave", "Taino", "Yukpa"),
                        pops_in_h3 = c("Andaman", "Australian", "Han", "French"),
                        pop_colors = c("blueviolet", "pink", 
                                       "darkgoldenrod1", "deepskyblue1", 
                                       "black"),
                        rm_legend = T){
  
  Dstat <- abba[(abba$H1 %in% pops_in_title |abba$H2 %in% pops_in_title) &
                  (abba$H3 %in% pops_in_h3),]
  df_y <- as.data.frame(pops_in_y)
  df_y$y_axis <- rownames(df_y)
  plot_list <- list()
  for(pop in pops_in_title){
    Dstat <- abba[(abba$H1 == pop |abba$H2 == pop) &
                    (abba$H3 %in% pops_in_h3),]
    Dstat <- switch_h2_h1(Dstat, pop)
    Dstat <- Dstat[!(Dstat$H1 %in% pops_in_h3),]
    
    # need this part for the gap in the plot
    fake_pop <- Dstat [1,]
    fake_pop[1, 1:8] <- 100
    fake_pop[1,9:12] <- c(pop, pop, 
                          pop, pop)
    Dstat <- rbind(Dstat, fake_pop)
    
    Dstat <- Dstat[Dstat$H1 %in% c(pop,pops_in_y),]
    Dstat$H1 <- factor(Dstat$H1, 
                        levels = rev(unique(c(pops_in_y, pop))))
    Dstat$H3 <- factor(Dstat$H3, 
                        levels = c(pops_in_h3, pop),
                        ordered = T)
    p <- ggplot(Dstat, aes(x = Z, y = H1, color = H3)) +
      
      coord_cartesian(xlim = c(-5,5)) +
      scale_color_manual(values = pop_colors,
                         breaks = pops_in_h3) +
      geom_vline(xintercept = -3.3, lty = "dashed", col = "gray")+
      geom_vline(xintercept = 3.3, lty = "dashed", col = "gray") +
      geom_vline(xintercept = 0, lty = "dashed", col = "gray") +
      geom_point(shape = 5, size = 3, stroke= 1) +
      labs(x = paste("D-stat(", pop, ", H2; H3, Yoruba)"), 
           title = pop) 
    
    if(rm_legend){
      p <- p + theme(legend.position = "none")
    }
    plot_list[[length(plot_list)+1]] <- p
  }
  return(plot_list)
}

prepare_set <- function(abba, pop_in_h2, pop_in_h3, pop_order){
  mySet <- abba[abba$H3 %in% pop_in_h3 ,]
  mySet <- switch_h2_h1(mySet, pop_in_h2 = pop_in_h2)
  mySet <- mySet[mySet$H2 == pop_in_h2,]
  mySet$H1 <- factor(mySet$H1, levels = pop_order, ordered = T)
  mySet <- mySet[order(mySet$H1),]
  rownames(mySet) <- seq(1, nrow(mySet))
  return(mySet)
}

# D/Z = std_err
get_stdErr <- function(D, Z){
  stdErr <- D/Z
  return(stdErr)
}

moreno_3c <- function(abba, title, myAlpha = 0.001,
                      xlim = c(-0.05, 0.05), yaxis = T,
                      color = "#d60036", cex = 1, lwd = 1,
                      pch = 5, col_to_plot = "JK.D", annotate = T){
  
  abba$stdErr <- qnorm(1 - myAlpha/2)*get_stdErr(D = abba[,col_to_plot],
                                                 Z = abba$Z)
  
  plot(x = abba[,col_to_plot], y = as.integer(rownames(abba)), 
       axes = F , xlab = "D-stat", ylab = NA,
       main = title, xlim = xlim,
       col = color, cex = cex, lwd = lwd,
       pch = pch)
  segments(y0 =  as.integer(rownames(abba)), 
           x0 = abba[,col_to_plot] - abba$stdErr,
           x1 = abba[, col_to_plot] + abba$stdErr,
           col = color
  )
  if(annotate){
  text(x = xlim[1]*0.8, 
       y = as.integer(rownames(abba)),
       labels = round(abba$Z,2))
  }
  grid(nx = NULL, ny = 0)
  if(yaxis){
    axis(2, at = rownames(abba), 
         labels = abba$H1, las = 2)
  }
  axis(1)
  abline(v = 0, lty = "dashed", col = "gray")
  
}

moreno_4a <- function(abba, pop, column_to_highlight, to_highlight,
                      color_to_highlight, default_color = "gray",
                      myAlpha = 0.001,
                      pch  = 5, cex = 1.2, lwd = 2, 
                      xaxis = T,
                      yaxis = T, plot = T, legend = F, order_h1 = FALSE,
                      thres_color = "black", x1 = F, x2 = F,
                      tick_color = F, groups_h1 = F, col_groups_h1 = F){
  
  # fix the colors of the dots
  abba$color <- default_color
  for(i in 1:length(to_highlight)){
    abba$color[abba[,column_to_highlight] == to_highlight[i]] <- color_to_highlight[i]
  }
  
  other_pops <- unique(abba[,column_to_highlight][!(abba[,column_to_highlight] %in% to_highlight)])
  
  abba[,column_to_highlight] <- factor(abba[,column_to_highlight], 
                                       levels = c(to_highlight, other_pops))
  
  # Order data frame to be plotted
  if(length(order_h1)>1){
    if(length(tick_color) > 1){
      abba[,col_groups_h1] <- factor(abba[,col_groups_h1], levels = groups_h1, ordered = T)
      
      abba <- abba[order(abba[, col_groups_h1], abba[,column_to_highlight],
                         decreasing = T),]
      abba$H1 <- factor(abba$H1, levels = unique(abba$H1), ordered = T)
      # Data frame with groups on y-axis and their color
      myTicks <- data.frame(col1 = groups_h1, col2 = as.character(tick_color))
      colnames(myTicks) <- c(col_groups_h1, "color_tick")
      abba <- join(abba, myTicks, by = col_groups_h1)
      # Data frame with populations on H1 (y-axis) and their group and tick color
      h1_tick <- unique(abba[,c("H1", col_groups_h1, "color_tick")])
      rownames(h1_tick) <- 1:nrow(h1_tick)
    }else{
      abba$H1 <- factor(abba$H1, levels = order_h1, ordered = T)
      abba <- abba[order(abba$H1,abba[,column_to_highlight],
                         decreasing = T),]
    }
  }else{
  abba <- abba[order(abba[,column_to_highlight], decreasing = T),]
  }
  rownames(abba) <- 1:nrow(abba)
  if(x1 & x2){
    xlim <- c(x1,x2)
  }else{
    xlim <- range(abba$Z)
  }
  if (plot){
    plot(x = abba$Z, y = abba$H1, col = abba$color, axes = F,
         xlab = "Z-score", ylab = NA,
         pch = pch, lwd = lwd, cex = cex, main = pop, xlim = xlim)
    if(yaxis){
      if(length(tick_color)>1){
        for(group in unique(h1_tick[,col_groups_h1])){
          myRows <- which(h1_tick[,col_groups_h1] == group)
          myCol <- unique(myTicks$color_tick[myTicks[, col_groups_h1] == group])
          axis(side=2, at = myRows,
               labels = h1_tick$H1[myRows], las = 2,
               col.axis = as.character(myCol))
        }
      }else{
        axis(side=2, at = 1:length(levels(abba$H1)),
             labels = levels(abba$H1), las = 2)
      }
      
    }
    if(xaxis){
      axis(1)
    }
    abline(v = c(-qnorm(1- myAlpha/2), qnorm(1- myAlpha/2)),
           lty = "dashed", col = thres_color)
  }
  if(legend){
    plot(1,1, type = "n", axes = F, xlab = NA, ylab = NA)
    legend("topleft",
           pch = pch, col = c(color_to_highlight, default_color),
           legend = c(to_highlight, "Other"), 
           lty = NA, bty = "n", title = "H3", cex = cex, lwd = lwd)
    
    legend("bottomleft",
           pch = NA, text.col = as.character(unique(h1_tick[,c(col_groups_h1, "color_tick")])$color_tick),
           legend = as.character(unique(h1_tick[,c(col_groups_h1)])),
           bty = "n")
  }
}

#plot_dstat <- function(abba, )