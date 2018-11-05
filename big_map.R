# Map for CIG retreat

library(ggplot2)
library(maptools )
library(stringi)
library(descr)
library(rgeos)
library(mapdata)

brMap <- readShapePoly("~/Projects/Botocudos/Files/Maps/estados_2010/estados_2010.shp")
minas <- brMap[brMap$nome=="Minas Gerais",]
espirito <- brMap[brMap$id == 8, ]
scatarina <- brMap[brMap$id == 24, ]
bahia <- brMap[brMap$id == 5,]

brMap$ESTADO<-
  stri_enc_toutf8(brMap$nome, T, T)
brMapDF <- as.data.frame(brMap) # para definir a region
brMap = gBuffer(brMap, width=0, byid=TRUE) #correct problem with Polygons - TopologyException
brMapDF

plot.map<- function(database,center,...){
  Obj <- map(database,...,plot=F)
  coord <- cbind(Obj[[1]],Obj[[2]])
  
  # split up the coordinates
  id <- rle(!is.na(coord[,1]))
  id <- matrix(c(1,cumsum(id$lengths)),ncol=2,byrow=T)
  if(dim(id)[1] > 1){
    polygons <- apply(id,1,function(i){coord[i[1]:i[2],]})
  }else{
    polygons <- list(coord)
  }

  
  # split up polygons that differ too much
  polygons <- lapply(polygons,function(x){
    if(!is.null(dim(x))){
      x[,1] <- x[,1] + center
      x[,1] <- ifelse(x[,1]>180,x[,1]-360,x[,1])
      if(sum(diff(x[,1])>300,na.rm=T) >0){
        id <- x[,1] < 0
        x <- rbind(x[id,],c(NA,NA),x[!id,])
      }
    }else{
      x <- x + center
      x <- ifelse(x>180,x-360,x)
      if(sum(diff(x)>300,na.rm=T) >0){
        id <- x < 0
        x <- rbind(x[id,],c(NA,NA),x[!id,])
      }
    }

    x
  })
  # reconstruct the object
  polygons <- do.call(rbind,polygons)
  Obj[[1]] <- polygons[,1]
  Obj[[2]] <- polygons[,2]
  
  map(Obj,...)
}

png("~/Projects/Botocudos/Plots/map_CIG.png",
    width = 20, height = 15, res = 350, units = "in")

plot.map("world", center=180, col="gray95",bg="cornflowerblue",
         fill=TRUE,ylim=c(-60,90),mar=c(0,0,0,0), border = F)

plot.map(database = "worldHires", regions = "Brazil", 
         center=180, col="black",
         fill=TRUE,ylim=c(-60,90),mar=c(0,0,0,0), add = T)

plot.map(minas, regions = "Minas Gerais",add=TRUE, 
     col=alpha("firebrick1", 0.7), 
     border=F, fill=TRUE,ylim=c(-60,90),mar=c(0,0,0,0),
     namefield = "nome", center = 180)

plot.map(espirito, add = T, 
     col = alpha("darkolivegreen1", 0.6), 
     border = F, fill=TRUE,ylim=c(-60,90),mar=c(0,0,0,0),
     namefield = "nome", center = 180)

plot.map(bahia, add = T, 
     col = alpha("gold", 0.6),
     border = F, fill=TRUE,ylim=c(-60,90),mar=c(0,0,0,0),
     namefield = "nome", center = 180)


dev.off()
