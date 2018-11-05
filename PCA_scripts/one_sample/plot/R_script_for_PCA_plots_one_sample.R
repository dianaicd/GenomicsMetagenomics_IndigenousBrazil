
setwd("/Users/vvilla/Dropbox/reporte_final/Tarahumara_remap/MOM6/PCA/")

sample="MOM6"

sample.pdf=paste(sample,"all_libraries.pdf",sep=".")

evec.in="MOM6_Complete_Final_Sequence.bwa.sort.rmdup.uniq.filtSites_MERGED.evec"
eval.in="MOM6_Complete_Final_Sequence.bwa.sort.rmdup.uniq.filtSites_MERGED.eval"
   
data <- read.table(evec.in)
eval <- read.table(eval.in)
#data <- data[nrow(data):1,]
## data must match the dimensions of data matrix (eg. (1,10) if the first 10 PCs are given) 
pca <- data[,seq(1,10)+1]

pct.varPC1 = paste("PC1 (",round(100*eval[1,]/sum(eval),2),"%)",sep="")
pct.varPC2 = paste("PC2 (",round(100*eval[2,]/sum(eval),2),"%)",sep="")
pct.varPC1
pct.varPC2

####
individual.name <- gsub(".*:","",data[,1])
sample.idx <- match(sample,gsub(".*:","",data[,1]))

#### population information:
popinfo <- read.csv("NatMex.txt",sep="\t",header=T)
## popinfo file can be longer than the actual number of samples ploted as long as the ID is matched by the index function
## but if there is an excess of samples with corresponding regions and colors, these will be ploted in the legend too
## therefore, ideal to have a popinfo file for each plotting session

regions <- as.character(unique(popinfo[,3]))
populations <- as.character(unique(popinfo[,5]))

regions.color <- c("black","#00008090", "#1E90FF90", "#0000FF90", "#87CEEB90", "#00F5F590")
#regions.color <- c("black","navyblue", "dodgerblue", "blue", "skyblue", "turquoise1", "darkgreen", "lightgreen")
regions.color <- c("black","#90ee9090", "#9932cc90", "#ff000090", "#ffa50090","#ffc0cb90") ## color of regions

names(regions.color) <- regions

print(cbind(regions,regions.color))

#### map to the pop infor.
nsamples <- length(individual.name)

idx <- match(individual.name, popinfo[,2])
pop.pc1.pc2 <- cbind(popinfo[idx,5],data[,2:3])
colors <- regions.color[as.character(popinfo[idx,3])]

##########################################################################################
## WITH SOLID CIRCLES
##########################################################################################

pdf(sample.pdf)
plot(pca[,1],pca[,2],col = colors, pch = 19,cex=1.3, xlab =pct.varPC1,ylab=pct.varPC2)
abline(v=0,h=0,lty=2,col="grey")
legend("topright", regions, ncol=1,col= regions.color,pch=19, cex=0.8,bty = "n")
points(pca[sample.idx,c(1,2)], cex = 1, pch = 19, col = "black")
text(pca[sample.idx,c(1,2)],sample,pos=4,offset=0.25,cex=0.6)

for (pop in populations){
  pop.idx <- grep(pop,pop.pc1.pc2[,1])
  pop.x <- mean (pop.pc1.pc2[pop.idx,]$V2)
  pop.y <- mean (pop.pc1.pc2[pop.idx,]$V3)
  text(pop.x,pop.y,pop,cex=0.4)
}


##########################################################################################
## WITH POP ID NAMES
##########################################################################################

#pdf(sample.ids.pdf)
plot(pca[,1],pca[,2],type="n",col = colors, pch = 19, xlab =pct.varPC1,ylab=pct.varPC2)
abline(v=0,h=0,lty=2,col="grey")
text(pca[,1],pca[,2], labels=popinfo[idx,1],col=colors,cex=0.5)
points(pca[sample.idx,c(1,2)], cex = 1, pch = 19, col = "black")
text(pca[sample.idx,c(1,2)],sample,pos=4,offset=0.25,cex=0.6)
dev.off()
