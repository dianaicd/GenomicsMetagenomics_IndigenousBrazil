#Rscript WorldHeatmap.R ../testinputfiles/Coords.txt 0 2 1.9 L mapita ALL
#Rscript WorldHeatmap.R ../testinputfiles/Coords.txt 0 2 1.9 R mapita -19.127125 34.310375 -34.995471 14.730652
#Rscript WorldHeatmap.R ../testinputfiles/Coords.txt 0 2 1.9 L mapita -19.127125 34.310375 -34.995471 14.730652
#Rscript WorldHeatmap.R ../testinputfiles/Coords.txt 0 2 1.9 L mapita ALL
#Rscript ~/Dropbox/PhD/WorldHeatmap/bin/WorldHeatmap_0_2.R f3\(YRI\;\ X\,\ 890\) 890_f3_coords.txt 0 3 1.9 R 890_f3_map ALL
#Rscript ~/Dropbox/PhD/WorldHeatmap/bin/WorldHeatmap_0_2.R f3\(YRI\;\ X\,\ 890\) 890_f3_coords.txt 1 3 1.9 R 890_f3_map ALL
#Rscript ~/Dropbox/PhD/WorldHeatmap/bin/WorldHeatmap_0_2.R f3\(YRI\;\ X\,\ 890\) 890_f3_coords.txt 1 3 1.9 R 890_f3_map -119.627263 -85 13.846990 33.402420 4 6
#Region+zoom: title infile center DecPl cex Side outprefix region regionName x1 x2 y1 y2
#Rscript ~/Dropbox/PhD/WorldHeatmap/bin/WorldHeatmap_0_3.R f3\(YRI\;\ X,\ Nr74\)\ Mean\ n.\ of\ sites\ =\ 271323 Nr74_f3_coords.txt 0 3 3 R Nr74_f3_map region Mexico -119.627263 -85 14 32 5 7

#side only works for zoom-ins. Can be L or anything else.


#InFile<-"../testinputfiles/Coords.txt"
#center<-0
#DecPl<-2
#OutFile<-"mapita"
#CexVal<-1.9
#SideForScale<-"L"

#x1<--19.127125
#x2<-34.310375
#y1<--34.995471
#y2<-14.730652

Args<-commandArgs(T)
Title<-Args[1]
Args<-Args[-1]
NumCategories<-as.numeric(Args[1])
Args<-Args[-1]
InFile<-Args[1]
center<-as.numeric(Args[2])
DecPl<-as.numeric(Args[3])
CexVal<-as.numeric(Args[4])
SideForScale<-Args[5]
OutFile<-Args[6]
zoomin<-F
if(Args[7]!="ALL"){
	if(Args[7]=="region"){
		library(maps)
		library(mapdata)
		RegionNames<-map("world", namesonly=TRUE, plot=FALSE)
		if(Args[8] %in% RegionNames){
			reg<-Args[8]
			heightVal<-as.numeric(Args[13])
			widthVal<-as.numeric(Args[14])
			
			x1<-as.numeric(Args[9])
			x2<-as.numeric(Args[10])
			y1<-as.numeric(Args[11])
			y2<-as.numeric(Args[12])

		}else{
			message("ERROR: Region no included in worldHires, check map(\"worldHires\", namesonly=TRUE, plot=FALSE)")
			q("no")
		}
	}else{
		zoomin<-T
		x1<-as.numeric(Args[7])
		x2<-as.numeric(Args[8])
		y1<-as.numeric(Args[9])
		y2<-as.numeric(Args[10])
		heightVal<-as.numeric(Args[11])
		widthVal<-as.numeric(Args[12])
	}
}

#center<-196
#DecPl<-2

library(maps)
library(mapdata)

#Awfully slow function from http://stackoverflow.com/questions/5353184/fixing-maps-library-data-for-pacific-centred-0-360-longitude-display
plot.map<- function(database,center,...){
    Obj <- map(database,...,plot=F)
    coord <- cbind(Obj[[1]],Obj[[2]])

    # split up the coordinates
    id <- rle(!is.na(coord[,1]))
    id <- matrix(c(1,cumsum(id$lengths)),ncol=2,byrow=T)
    polygons <- apply(id,1,function(i){coord[i[1]:i[2],]})

    # split up polygons that differ too much
    polygons <- lapply(polygons,function(x){
        x[,1] <- x[,1] + center
        x[,1] <- ifelse(x[,1]>180,x[,1]-360,x[,1])
        if(sum(diff(x[,1])>300,na.rm=T) >0){
          id <- x[,1] < 0
          x <- rbind(x[id,],c(NA,NA),x[!id,])
       }
       x
    })
    # reconstruct the object
    polygons <- do.call(rbind,polygons)
    Obj[[1]] <- polygons[,1]
    Obj[[2]] <- polygons[,2]

    map(Obj,...)
}

transcol<-Vectorize(function(x, a){
	col<-col2rgb(x)
	newcol<-rgb(col[1]/255, col[2]/255, col[3]/255, a)
	return(newcol)
})


a<-read.table(InFile)
if(Args[7]=="region" | zoomin==T){
	a<-a[!(a[,1]<y1 | a[,1]>y2 | a[,2]<x1 | a[,2]>x2),]
	#save(a, y1, y2, x1, x2, file="sisisi.rd")
}
a<-a[sort.list(a[,3], decreasing=T),]
b<-a[,2]
b<-b+center
b<-ifelse(b>180,b-360,b)


cols<-transcol(colorRampPalette(c("blue4", "darkblue", rev(rainbow(15)[-c(1,5,6,7,12,13,14,15)]), "red"))(NumCategories), .75)

z=matrix(1:NumCategories,nrow=1)
x=1
y=seq(min(a[,3], na.rm=T),max(a[,3], na.rm=T),len=NumCategories)

my.colors<-transcol(colorRampPalette(c("blue4", "darkblue", rev(rainbow(15)[-c(1,5,6,7,12,13,14,15)]), "red"))(NumCategories), .75)

levs<-seq(min(a[,3], na.rm=T), max(a[,3], na.rm=T), (max(a[,3], na.rm=T)-min(a[,3], na.rm=T))/(NumCategories-1))
newcols<-NULL
for(i in a[,3]){
	newcols<-c(newcols, cols[which.min(abs(levs-i))])
}




if(center==0){
message("Plotting Atlantic-centered map... ")

if(zoomin){

message("Plotting zoom-in... ")

#pdf(paste(OutFile, ".pdf", sep=""), height=abs(y2-y1)*3/180*2, width=abs(x2-x1)*7/360*1.7)
#pdf(paste(OutFile, ".pdf", sep=""), height=3.5, width=6)
pdf(paste(OutFile, ".pdf", sep=""), height=heightVal, width=widthVal)
if(SideForScale=="L"){
layout(matrix(c(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1), 1, 22, byrow = TRUE))
}else{
layout(matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2), 1, 22, byrow = TRUE))
}
par(fg="white")
#plot.map("worldHires", col="white", fill=T, bg="gray85", mar=c(0, 0, 0, 0),ylim=c(y1,y2), xlim=c(x1,x2), center=center)
map("world", col="white", fill=T, bg="gray85", mar=c(0, 0, 0, 0),ylim=c(y1,y2), xlim=c(x1,x2))
#points(a[,2], a[,1], col="black", pch=16, cex=1.6)
points(a[,2], a[,1], col=newcols, pch=16, cex=CexVal)
#points(b, a[,1], col=cols, pch=16, cex=1.5)
title(paste(sep="", "\n", Title))
par(mar=c(1,0,1,0))
par(fg="darkblue")
image(x,y,z,col=my.colors,axes=FALSE,xlab="",ylab="")
if(SideForScale=="L"){
axis(4, at=seq(min(a[,3], na.rm=T), max(a[,3], na.rm=T), (max(a[,3], na.rm=T)-min(a[,3], na.rm=T))/5), labels=round(seq(min(a[,3], na.rm=T), max(a[,3], na.rm=T), (max(a[,3], na.rm=T)-min(a[,3], na.rm=T))/5), DecPl), las=2)
}else{
axis(2, at=seq(min(a[,3], na.rm=T), max(a[,3], na.rm=T), (max(a[,3], na.rm=T)-min(a[,3], na.rm=T))/5), labels=round(seq(min(a[,3], na.rm=T), max(a[,3], na.rm=T), (max(a[,3], na.rm=T)-min(a[,3], na.rm=T))/5), DecPl), las=2)
}
dev.off()


}else if(Args[7]=="region"){

message("Plotting region zoom-in... ")

pdf(paste(OutFile, ".pdf", sep=""), height=heightVal, width=widthVal)
layout(matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2), 1, 22, byrow = TRUE))
par(fg="white")
map("world", reg, col="white", fill=T, bg="gray85", ylim=c(y1,y2), xlim=c(x1,x2), mar=c(0,0,1,0))
title(Title)
#points(a[,2], a[,1], col="black", pch=16, cex=1.6)
points(a[,2], a[,1], col=newcols, pch=16, cex=CexVal)
#points(b, a[,1], col=cols, pch=16, cex=1.5)
par(mar=c(5,0,5,0))
par(fg="darkblue")
image(x,y,z,col=my.colors,axes=FALSE,xlab="",ylab="")
axis(2, at=seq(min(a[,3]), max(a[,3]), (max(a[,3])-min(a[,3]))/5), labels=round(seq(min(a[,3]), max(a[,3]), (max(a[,3])-min(a[,3]))/5), DecPl), las=2)
dev.off()

}else{


pdf(paste(OutFile, ".pdf", sep=""), height=3)
layout(matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2), 1, 22, byrow = TRUE))
par(fg="white")
#plot.map("worldHires", col="white", fill=T, bg="gray85", mar=c(0, 0, 0, 0),ylim=c(-60,90), center=center)
map("world", col="white", fill=T, bg="gray85", mar=c(0, 0, 0, 0),ylim=c(-60,90))
title(Title)
#points(a[,2], a[,1], col="black", pch=1, cex=CexVal)
points(a[,2], a[,1], col=newcols, pch=16, cex=CexVal)
#points(b, a[,1], col=cols, pch=16, cex=1.5)
par(mar=c(5,0,5,0))
par(fg="darkblue")
image(x,y,z,col=my.colors,axes=FALSE,xlab="",ylab="")
axis(2, at=seq(min(a[,3]), max(a[,3]), (max(a[,3])-min(a[,3]))/5), labels=round(seq(min(a[,3]), max(a[,3]), (max(a[,3])-min(a[,3]))/5), DecPl), las=2)
dev.off()

}

}else{
message("Plotting Pacific-centered map... ")

pdf(paste(OutFile, ".pdf", sep=""), height=3)
layout(matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2), 1, 22, byrow = TRUE))
par(fg="white")
plot.map("world", col="white", fill=T, bg="gray85", mar=c(0, 0, 0, 0),ylim=c(-60,90), center=center)
title(paste(sep="", "\n", Title))
#map("worldHires", col="white", fill=T, bg="gray85", mar=c(0, 0, 0, 0),ylim=c(-60,90))
#points(a[,2], a[,1], col="black", pch=16, cex=1.6)
#points(a[,2], a[,1], col=cols, pch=16, cex=1.5)
points(b, a[,1], col=newcols, pch=16, cex=CexVal)
par(mar=c(5,0,5,0))
par(fg="darkblue")
image(x,y,z,col=my.colors,axes=FALSE,xlab="",ylab="")
axis(2, at=seq(min(a[,3]), max(a[,3]), (max(a[,3])-min(a[,3]))/5), labels=round(seq(min(a[,3]), max(a[,3]), (max(a[,3])-min(a[,3]))/5), DecPl), las=2)
dev.off()

}



