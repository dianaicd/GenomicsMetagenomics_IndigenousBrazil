mds_plot <- function(ind, d, ndim = 2, colors = F, cx = 1, cy =  2, label = F, title = F){
  panel <- read.table("~/Projects/Botocudos/Files/PCA/2017_09_15/panel_botos.txt",
                      header = T)
  
  path <- paste(d, ind, sep = "")
  distance <- read.table(paste(path, ".dist", sep = ""))
  ids <- read.table(paste(path, ".dist.id", sep = ""))
  
  result <- cmdscale(distance, k = ndim, x.ret = T)
  coordinates <- as.data.frame(result$points)
  coordinates <- cbind(coordinates, ids)
  colnames(coordinates) <- c("x", "y", "region", "indID")
  coordinates$region<- NULL
  coordinates$indID <- sub("-", "_", coordinates$indID)
  coordinates <- join(coordinates, panel, by = "indID")
  l <- levels(coordinates$population)
  index <- grepl("Boto", l)
  levels(coordinates$population) <- c("Botocudos", l[!index])
#  coordinates$region[grep("MA23", coordinates$population)] <- "Botocudo"
  if(!is.character(colors)){
    colors <- funky(length(unique(coordinates$population))-1)
    colors <- c("black", colors)
  }
  if(!is.character(title)){
    title <- ind
  }
  
  if(label){
    p <- ggplot(coordinates[rev(rownames(coordinates)),], 
                +aes(x = x, y = y, 
                     color = population, label = population)) +
      geom_text() +
      scale_color_manual(values = colors, name = "Region") +
      labs(x = paste("Coordinate", cx), y = paste("Coordinate", cy),
           title = title)
  }else{

    p <- ggplot(coordinates[rev(rownames(coordinates)),], 
                aes(x = x, y = y, color = population)) +
      geom_point(size=4) +
      scale_color_manual(values = colors, name = "Region") +
      labs(x = paste("Coordinate", cx), y = paste("Coordinate", cy),
           title = title)
  }
  
  return(p)
}
