order_levels <- function(coord, colorby, divisions, division_type){
  coord[,colorby] <- droplevels(coord[,colorby])
  old_levels <- levels(coord[,colorby])
  ordered_levels <- c(old_levels[!old_levels %in% divisions], divisions)
  ordered_levels <- ordered_levels[ordered_levels %in% old_levels]
  index <- unlist(
    sapply(ordered_levels, 
           function(x) 
             seq(1, dim(coord)[1])[grepl(x, coord[,division_type])]))
  coord <- coord[index, ]

  coord[,division_type] <- factor(coord[,division_type], 
                                  levels =  ordered_levels)
  
  return(coord)
  # l <- levels(coord[,colorby])
  # new_levels <- c()
  # for(i in levels(coord[,colorby])){
  #   index <- grepl(i, coord[,colorby])
  #   new_levels <- c(new_levels, as.character(coord[index,colorby]))
  # }
  # index <- coord[,colorby] %in% new_levels
  # new_levels <- c(new_levels, as.character(coord[!index, colorby]))
  # #  index <- grep("Boto", l)
  # coord[,colorby] <- factor(coord[,colorby], levels =  unique(new_levels))
  # index <- levels(coord[,colorby]) %in% coord[, colorby]
}

substitution <- function(panel, coord, substitute){
  # A silly solution for the different ways in which people
  # labeled individuals in the panels
  panel$indID <- gsub("[[:punct:]]", "_", panel$indID)
  coord$Individual <- gsub("[[:punct:]]", "_", coord$Individual)
  coord$pop_label <- coord$population
  if(substitute){
    for(i in panel$indID){
      ex <- paste("[[:alnum:]]*_?", i, "_?[[:alnum:]]*", sep = "")
      coord$Individual <- gsub(ex, i, coord$Individual)
    }
    coord$Individual <- gsub("Coastal_", "", coord$Individual)
    coord$Individual <- gsub("USAmerindian_", "", coord$Individual)
    coord$Individual <- gsub("CanAmerindian_", "", coord$Individual)
  }
  
  coord$indID <- coord$Individual
  coord <- join(coord, panel, by = "indID")
  
  return(coord)
  
}

plot_colors <- function(colors, coord, colorby){
  index <- levels(coord[,colorby]) %in% coord[, colorby]
  
  if(!is.character(colors)){
    colors <- funky(sum(index)-1)
    #lab_col <- cbind(c(0,0,0), lab_col  - min(lab_col))
    colors <- c(colors, "black")
    lab_col <- col2rgb(colors)
    names(colors) <- levels(coord[,colorby])
  }else{
    lab_col <- col2rgb(colors)
  }
  
  # This part is to make the colors of the labels a bit darker
  # than the colors for the dots
  lab_col[lab_col > 0] <- lab_col[lab_col > 0]  - min(lab_col[lab_col > 0])
  lab_col <- rgb(t(lab_col)/255)
  
  lab_col <- data.frame(Color = lab_col, 
                        colorby = levels(coord[, colorby])[index])
  lab_col <- lab_col[order(lab_col$colorby),]
  return(list(colors, lab_col))
}



plot_bammds <- function(d, ind, colors = F, label = F,
                        title=F, cx = 1, cy = 2, 
                        divisions = c("Botocudos", "America", "Oceania"),
                        switch = F,
                        colorby = "population",
                        panel_path = "~/Projects/Botocudos/Files/PCA/2017_09_15/panel_botos.txt",
                        substitute=T,
                        addlabel=F,
                        division_type=colorby){

  panel <- read.table(panel_path,
                      header = T, sep = "\t")
  
  coord_path <- paste(d, ind, "_d", cx, "_d", cy, ".csv", sep = "")
  coord <- read.csv(coord_path, header = T)
  coord <- substitution(panel, coord, substitute)
  coord <- order_levels(coord, colorby, divisions, division_type)
  tmp_color <- plot_colors(colors, coord, colorby)
  colors <- tmp_color[[1]]
  lab_col <- tmp_color[[2]]
    # Plot title
  if(!is.character(title)){
    title <- ind
  }
  
  #coord <- coord[order(coord$population, decreasing = T),]
  # annotation
  annotation <- data.frame()
  for(pop in levels(coord[, colorby])){
    L <- length(coord[, colorby])
    index <- grepl(pop, coord[, colorby])
    # Make sure that we are sampling
    index <- sample(x = rep(seq(1, L)[index], 2), 1)
    annotation <- rbind(annotation, coord[index, ])
  }
  annotation <- annotation[order(as.character(annotation[,colorby])),]
  
  if(switch){
    tmp <- cx
    cx <- cy
    cy <- tmp
    rm(tmp)
    colnames(coord) <- sub("mds_x", "tmp", colnames(coord))
    colnames(coord) <- sub("mds_y", "mds_x", colnames(coord))
    colnames(coord) <- sub("tmp", "mds_y", colnames(coord))
    colnames(annotation) <- colnames(coord)
  }
  
  if(label){
    p <- ggplot(coord, 
                aes(x = mds_x, y = mds_y, 
                    color = coord[,colorby], label = coord[,colorby])) +
      geom_text() +
      scale_color_manual(values = colors, name = "Region") +
      labs(x = paste("Coordinate", cx), y = paste("Coordinate", cy),
           title = title) 
    L <- get_legend(p)
    p <- p + theme(legend.position = "none")
    plots <- list(p, L)
  }else{
    p <- ggplot(coord, 
                aes(x = mds_x, y = mds_y, color = coord[,colorby])) +
      geom_point(size=6, alpha=0.7) +
      scale_color_manual(values = colors, name = "Region") +
      labs(x = paste("Dimension", cx), y = paste("Dimension", cy),
           title = title)
    
    if(addlabel){
      p <- p + geom_text(
        data = annotation, 
        aes(x = mds_x, y = mds_y), color = lab_col$Color, 
        label = annotation[, colorby], 
        stat = "identity")
    }

    L <- get_legend(p)
    p <- p + theme(legend.position = "none")
    plots <- list(p, L)
  }
  
   return(plots)
}
