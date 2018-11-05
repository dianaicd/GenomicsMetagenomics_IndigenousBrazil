build_sex <- function(d="~/Projects/Botocudos/Files/Sex_determination/2017_09_15/", boto){
  setwd(d)
  sex <- data.frame()
  files <- list.files(pattern = "*table")
  for(f in files){
    tmp <- read.table(f, header = T, sep = "\t")
    ind <- gsub("\\..*", "", f)
    tmp$Library <- ind
    sex <- rbind(sex, tmp)
  }
  sex$Library <- factor(sapply(sex$Library, function(x) ma2mn(x)),
                        levels = order_endo(boto=boto)$Library)
                          # c("MN00010", "MN00013", "MN00016", "MN00019", "MN00021", "MN00022", "MN00023", "MN0003", 
                                                       #"MN00039", "MN00045", "MN00056", "MN00064", "MN00066", "MN00067", "MN00068", "MN00069",
                                                       #"MN0009",  "MN00118", "MN00119", "MN00316", "MN00346", "MN01701", "MN1943" ))
  return(sex)
}


ggsex <- function(d = "~/Projects/Botocudos/Files/Sex_determination/2017_09_15/", boto){
  sex <- build_sex(d=d, boto = boto)
  colors <- c("palevioletred4", "royalblue4", "royalblue1", "palevioletred2")
  names(colors) <- levels(sex$Assignment)
  
  ggplot(sex, aes(x = Library, y = R_y,
                  ymin = R_y - 2*SE, ymax = R_y + 2*SE,
                  color = Assignment)) +
    coord_flip(xlim = c(-1,dim(sex)[1]+2)) +
    geom_point() +
    geom_errorbar() +
    labs(y = expression(R[y]), x = NULL, title = "Sex determination") +
    geom_hline(color = "gray", size = 0.8, yintercept = 0.042, lty = "dashed") +
    scale_color_manual(values = colors) +
    geom_rect(aes(ymin = 0, ymax = 0.016, xmin = 0, xmax = 24), color = "snow2",
              fill = "snow2") +
    geom_rect(aes(ymin = 0.075, ymax = 0.11, xmin = 0, xmax = 24), color = "snow2",
              fill = "snow2") +
    geom_point() +
    geom_errorbar() +
    theme(legend.position = "right") +
    annotate("text", y = 0.01, x = dim(sex)[1]+1, label = "Female") +
    annotate("text", y = 0.09, x = dim(sex)[1]+1, label = "Male")
    
}
