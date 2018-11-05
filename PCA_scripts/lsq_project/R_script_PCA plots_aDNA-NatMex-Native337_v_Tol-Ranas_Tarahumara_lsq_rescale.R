
### For PCA plot of F9_MOM6_lsq
#setwd("/Users/vvilla/Desktop/plost_Without_LAC-SER/lsq/Tarahumara_Tol-Ranas/")

setwd("~/Dropbox (UNAM-LIIGH)/LIIGH/STUDENTS/PhDStudents/Viridiana/ProjectoTarahumara/PCA_lsq_Tarahumara_Tol-Ranas")

sample="Tol-Ranas_Tarahumara_lsq_rescale_5"

sample.svg=paste(sample,"3b.svg",sep=".")

evec.in="Tol-Ranas_Tarahumara_NatAmMex_without_LAC-SER_5pop.lsq.evec"
eval.in="Tol-Ranas_Tarahumara_NatAmMex_without_LAC-SER_5pop.lsq.eval"

data <- read.table(evec.in, col.names= c("sample", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "Pop"))
eval <- read.table(eval.in)
#data <- data[nrow(data):1,]
## data must match the dimensions of data matrix (eg. (1,10) if the first 10 PCs are given) 
## seq (Sequence generation. Generate regular sequences)
### usage: seq(from = 1, to = 1, by = ((to - from)/(length.out - 1)),
#####		    length.out = NULL, along.with = NULL, ...)
pca <- data[,seq(1,10)+1]

### paste (concatenate strings. concatenate vectors after converting to character)
### round (Rounding numbers. round rounds the values in its first argument to the specified number of decimal places.
#### usage: round(x, digits = 0))
pct.varPC1 = paste("PC1 (",round(100*eval[1,]/sum(eval),2),"%)",sep="")
pct.varPC2 = paste("PC2 (",round(100*eval[2,]/sum(eval),2),"%)",sep="")
pct.varPC1
pct.varPC2

#### gsub (Pattern matching and replacement. usage: gsub(pattern, replacement, x, ignore.case = FALSE, perl = FALSE,
###     fixed = FALSE, useBytes = FALSE))
individual.name <- gsub(".*:","",data[,1])

### match (value matching. match returns a vector of the positions of first matches 
#####of its first argument in its second. usage: match(x, table, nomatch = NA_integer_, incomparables = NULL))
sample.idx <- match(sample,gsub(".*:","",data[,1])) #### ok indica que se encuentra en la fila 1


#### population information:
popinfo <- read.csv("Ancient-Tol-Ranas_Tarahumara_NatAM_MEX_demodata_without_LAC-SER_2.txt",sep="\t",header=T)
## popinfo file can be longer than the actual number of samples ploted as long as the ID is matched by the index function
## but if there is an excess of samples with corresponding regions and colors, these will be ploted in the legend too
## therefore, ideal to have a popinfo file for each plotting session


regions <- as.character(unique(popinfo[,3]))
populations <- as.character(unique(popinfo[,5]))

regions.color <- c("darkred","deeppink","black","navyblue", "deepskyblue","#90ee9090", "#9932cc90", "#ff000090", "#ffa50090","#ffc0cb90")
####Cambiar colores para que sea mas clara la diferencia entre las poblaciones

names(regions.color) <- regions ### ok


#### cbind ( combine R objects by rows or columns.Take a sequence of vector, matrix or data frames arguments and combine 
###by columns or rows, respectively. These are generic functions with methods for other R classes. 
#####usage: cbind(..., deparse.level = 1) ) (Arguments: ... vector or matrices)
print(cbind(regions,regions.color))

#### map to the pop infor.
### length (length of an object. get or set the length of vectors including lists and factors, and of any other R object 
###for which a method has been defined. usage: length(x))
nsamples <- length(individual.name)

idx <- match(individual.name, popinfo[,2])
pop.pc1.pc2 <- cbind(popinfo[idx,5],data[,2:3])
colors <- regions.color[as.character(popinfo[idx,3])]

shape <- c(8,13,11,9, 14, 19, 19, 19, 19,19, 19, 19, 19, 19,19,19,19)
names(shape) <- populations
pch.shape <- shape[as.character(popinfo[idx,5])]

shape2 <- c(8,13,11,9, 14, 19, 19, 19, 19,19)
#names(shape2) <- regions
#pch.shape2 <- shape2[as.character(popinfo[idx,3])]



##########################################################################################
## WITH SOLID CIRCLES
##########################################################################################

### pdf (PDF graphics device, pdf graphics device driver for producng PDF grahics)
svg(sample.svg)
## pch: ploting character , cex: character expansion, a numerical vector.

pdf("Plot_4_cancun_maria.pdf")
plot(pca[,1],pca[,2],col = colors, pch = pch.shape,cex=1.2, xlab =pct.varPC1,ylab=pct.varPC2)
## lty: line type, 2= dashed
abline(v=0,h=0,lty=2,col="grey")
### ncol: the number of columns in which to set the legend items, the default is 1, a vertical legend.
####bty: a character string which determined the type of box which is drawn about plots. 
legend("bottomright", regions, ncol=1,col= regions.color,pch= shape2, cex=0.8,bty = "n")
#### points (add points to a plot. points is a generic function to draw a sequence of points at the specified coordinates.
### The specified characters are plotted, centered at the coordinates)
points(pca[sample.idx,c(1,2)], cex = 1, pch = 19, col = "black")
### text. add text to a plot, text drawss the strings given in the vector labels at the coordinates given by x and y. 
### "y" may be missing since xy coords is used for construction of the coordinates.
### pos: a position specifier for the text. if specified this overrides any adj value given. Values of 1, 2, 3 and 4, 
### respectively indicate positions below, to the left of, above and to the right of the specified coordinates. 
### offset: when pos is specified, this value gives the offset of the label from the specified coordinate in fractions 
### of a character width.
dev.off()
#text(pca[sample.idx,c(1,2)],sample,pos=4,offset=0.25,cex=0.6)

#for (pop in populations){
  #pop.idx <- grep(pop,pop.pc1.pc2[,1])
  #pop.x <- mean (pop.pc1.pc2[pop.idx,]$V2)
  #pop.y <- mean (pop.pc1.pc2[pop.idx,]$V3)
  #text(pop.x,pop.y,pop,cex=0.4)
#}


##########################################################################################
## WITH POP ID NAMES
##########################################################################################

#pdf(sample.ids.pdf)
#plot(pca[,1],pca[,2],type="n",col = colors, pch = 19, xlab =pct.varPC1,ylab=pct.varPC2, sub= "lsqproject")
#abline(v=0,h=0,lty=2,col="grey")
#text(pca[,1],pca[,2], labels=popinfo[idx,5],col=colors,cex=0.8)
#points(pca[sample.idx,c(1,2)], cex = 1, pch = 17, col = "black")
#text(pca[sample.idx,c(1,2)],sample,pos=3,offset=0.2,cex=0.8)
#legend("bottomright", regions, ncol=1,col= regions.color,pch=19, cex=0.8,bty = "n")
#dev.off()













