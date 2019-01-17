library(mapdata)

# limits for latitude
y1 <- 30
y2 <- 45
#limits for longitude
x1 <- 18
x2 <- 40

# Sea color
sea <- "powderblue"
# Land color
land <- "white"
border <- "gray"

# Name of file and resolution, height, width.
png("~/Desktop/greece.png", width = 5, height = 5, res = 250, units = "in")

# Plot Greece
map(database = "worldHires", regions = "Greece", 
         col=land,
         fill=TRUE,ylim=c(y1,y2),mar=c(0.1,0.1,0.1,0.1), bg = sea,
         border = border, xlim = c(x1, x2))
# Add other countries
map(database = "worldHires", regions = "Turkey", 
          col=land,
         fill=TRUE,ylim=c(y1,y2),
         border = border, add = T)
map(database = "worldHires", regions = "Bulgaria", 
          col="white",
         fill=TRUE,ylim=c(y1,y2),mar=c(0,0,0,0),
         border = border, add = T)
map(database = "worldHires", regions = "Albania", 
         col="white",
         fill=TRUE,ylim=c(y1,y2), bg = "powderblue",
         border = border, add = T)
map(database = "worldHires", regions = "Yugoslavia", 
         col="white",
         fill=TRUE,ylim=c(y1,y2), bg = "powderblue",
         border = border, add = T)
map(database = "worldHires", regions = "Italy", 
    col=land,
    fill=TRUE,ylim=c(y1,y2), bg = "powderblue",
    border = border, add = T)

dev.off()

