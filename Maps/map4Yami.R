
library(ggplot2)
library(maptools )
library(stringi)
library(descr)
library(rgeos)
library(mapdata)
library(rgdal)

brMap <- readOGR("~/Projects/Botocudos/Files/Maps/estados_2010/estados_2010.shp")
minas <- brMap[brMap$nome=="Minas Gerais",]
espirito <- brMap[brMap$id == 8, ]
scatarina <- brMap[brMap$id == 24, ]
bahia <- brMap[brMap$id == 5,]

brMap$ESTADO<-
  stri_enc_toutf8(brMap$nome, T, T)
brMapDF <- as.data.frame(brMap) # para definir a region
brMap = gBuffer(brMap, width=0, byid=TRUE) #correct problem with Polygons - TopologyException
brMapDF


library(mapdata)

# png("~/Projects/Botocudos/Plots/map_colombia.png",
#     width = 5, height = 5, res = 200, units = "in")
map("worldHires","Brazil", col="black", fill=TRUE,
    bg = "cornflowerblue",
    xlim = c(-95, -25), ylim = c(-35, 20), mar = c(0,0,0,0)
)

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
#dev.off()
