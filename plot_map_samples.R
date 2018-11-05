library(ggplot2)
#library(ggvis)
#library(ggmap)
library(maptools )
library(stringi)
library(descr)
#library(rgdal)
library(rgeos)
#library(mapplots)
#library(maps)
library(mapdata)
brMap <- readShapePoly("~/Projects/Botocudos/Files/Maps/estados_2010/estados_2010.shp")
minas <- brMap[brMap$nome=="Minas Gerais",]
espirito <- brMap[brMap$id == 8, ]
scatarina <- brMap[brMap$id == 24, ]
bahia <- brMap[brMap$id == 5,]
#brMap
brMap$ESTADO<-
  stri_enc_toutf8(brMap$nome, T, T)
brMapDF <- as.data.frame(brMap) # para definir a region
brMap = gBuffer(brMap, width=0, byid=TRUE) #correct problem with Polygons - TopologyException
brMapDF

pdf("../Plots/Brazil_map.pdf")
png("~/Projects/Botocudos/Plots/Brazil_map_points.png", res = 200, width = 4, height = 4, units = "in")
map("worldHires","Brazil", col="black", fill=TRUE,
    bg = "azure",
    xlim = c(-75, -35), ylim = c(-35, 10), mar = c(0,0,0,0)
    #xlim = c(-52, -38), ylim = c(-25, -13)
    )
plot(minas, add=TRUE, 
     col=alpha("firebrick1", 0.7), 
     border=F)
plot(espirito, add = T, 
     col = alpha("darkolivegreen1", 0.6), 
     border = F)
plot(scatarina, add = T, 
     col = alpha("hotpink1", 0.6),
     border = F)
plot(bahia, add = T, 
     col = alpha("gold", 0.6),
     border = F)
points(x = c(-45, -41, -50, -42), 
       y = c(-19, -20, -27, -12), 
       col = c(alpha("white",0.7), 
               alpha("white", 0.6), 
               alpha("white", 0.6),
               alpha("white", 0.6)),
       cex = c(10/8,14/8, 5/8,5/8), pch = 16)
dev.off()

map("worldHires","venezuela", col="gray95",fill=TRUE, add=TRUE, border = F)  #add the adjacent parts of the US; can't forget my homeland
map("worldHires","chile", col="gray95",fill=TRUE, add=TRUE, border = F)
map("worldHires","colombia", col="gray95",fill=TRUE, add=TRUE, border = F)
map("worldHires","guyana", col="gray95", fill=TRUE, add=TRUE, border = F)
map("worldHires","surinam", col="gray95", fill=TRUE, add=TRUE, border = F)
map("worldHires","peru", col="gray95", fill=TRUE, add=TRUE, border = F)
map("worldHires","bolivia", col="gray95", fill=TRUE, add=TRUE, border = F)
map("worldHires","paraguay", col="gray95", fill=TRUE, add=TRUE, border = F)
map("worldHires","argentina", col="gray95", fill=TRUE, add=TRUE, border = F)
map("worldHires","uruguay", col="gray95", fill=TRUE, add=TRUE, border = F)
map("worldHires","french guiana", col="gray95", fill=TRUE, add=TRUE, border = F)
dev.off()

add.pie(z=c(1), x=-44, y=-18, radius=sqrt(10), col=alpha("gold", 0.6), labels="")

#z indicates the portions of the pie charts filled by each given type, 
#x & y are coordinates for the point, and radius is to designate the size of the circle for the pie chart


#to plot all points, I run a loop to run through my data one point at a time and make each pie chart; there are most likely more efficient methods

#######
# Zoom
map("worldHires","Brazil", col="black", fill=TRUE,
    bg = "cornflowerblue",
    #xlim = c(-75, -35), ylim = c(-35, 10)
    xlim = c(-57, -38), ylim = c(-35, -13)
)
plot(minas, add=TRUE, 
     col=alpha("gray90", 0.1), 
     border="firebrick1")
plot(espirito, add = T, 
     col = alpha("gray90", 0.1), 
     border = "darkolivegreen1")
plot(scatarina, add = T, 
     col = alpha("gray90", 0.1),
     border = "hotpink1")
points(x = c(-45, -41, -50), 
       y = c(-19, -20, -27), 
       col = c(alpha("firebrick1",0.7), 
               alpha("darkolivegreen1", 0.6), 
               alpha("hotpink1", 0.6)),
       cex = c(4,3,1), pch = 16)

add.pie(z=c(1), x=-44, y=-18, radius=sqrt(3), 
        col=alpha("gold", 0.9),
        labels="")


#### Nicer map
png("~/Projects/Botocudos/Plots/map_CIG.png",
    width = 7, height = 7, res = 350, units = "in")

map("worldHires","Brazil", col="black", fill=TRUE,
    bg = "cornflowerblue",
    xlim = c(-95, -25), ylim = c(-35, 20), mar = c(0,0,0,0)
    #xlim = c(-52, -38), ylim = c(-25, -13)
)

plot(minas, add=TRUE, 
     col=alpha("firebrick1", 0.7), 
     border=F)
plot(espirito, add = T, 
     col = alpha("darkolivegreen1", 0.6), 
     border = F)
plot(scatarina, add = T, 
     col = alpha("hotpink1", 0.6),
     border = F)
plot(bahia, add = T, 
     col = alpha("gold", 0.6),
     border = F)
points(x = c(-45, -41, -50, -42), 
       y = c(-19, -20, -27, -12), 
       col = c(alpha("white",0.7), 
               alpha("white", 0.6), 
               alpha("white", 0.6),
               alpha("white", 0.6)),
       cex = c(10/8,14/8, 5/8,5/8), pch = 16)
map("worldHires","venezuela", col="gray95",fill=TRUE, add=TRUE, border = F)  #add the adjacent parts of the US; can't forget my homeland
map("worldHires","chile", col="gray95",fill=TRUE, add=TRUE, border = F)
map("worldHires","colombia", col="gray95",fill=TRUE, add=TRUE, border = F)
map("worldHires","guyana", col="gray95", fill=TRUE, add=TRUE, border = F)
map("worldHires","surinam", col="gray95", fill=TRUE, add=TRUE, border = F)
map("worldHires","peru", col="gray95", fill=TRUE, add=TRUE, border = F)
map("worldHires","bolivia", col="gray95", fill=TRUE, add=TRUE, border = F)
map("worldHires","paraguay", col="gray95", fill=TRUE, add=TRUE, border = F)
map("worldHires","argentina", col="gray95", fill=TRUE, add=TRUE, border = F)
map("worldHires","uruguay", col="gray95", fill=TRUE, add=TRUE, border = F)
map("worldHires","french guiana", col="gray95", fill=TRUE, add=TRUE, border = F)
map("worldHires","ecuador", col="gray95", fill=TRUE, add=TRUE, border = F)
map("worldHires","panama", col="gray95", fill=TRUE, add=TRUE, border = F)
map("worldHires","costa rica", col="gray95", fill=TRUE, add=TRUE, border = F)
map("worldHires","nicaragua", col="gray95", fill=TRUE, add=TRUE, border = F)
map("worldHires","honduras", col="gray95", fill=TRUE, add=TRUE, border = F)
dev.off()
