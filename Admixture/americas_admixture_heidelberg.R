library(plyr)
setwd("~/Projects/Botocudos/Files/ADMIX/2019_03_17/Worldwide/")
panel <- read.table("~/Projects/Botocudos/Files/ADMIX/2019_03_17/Maanasa_0.01_minind_americas_Botocudos_fromVictor.panel", header = T)
sufix <- "Subset_Botocudos_fromVictor_k"

k <- "2_2"
name <- paste(sufix, k, ".qopt", sep = "")
admix <- t(as.matrix(read.table(name)))
nInd <- dim(admix)[2]



# panel$Region <- factor(panel$Region, levels = c("Africa", "Europe", "Caucasus",
#                                                 "NearEast", "SouthAsia",
#                                                 "SoutheastAsia", "EastAsia", "Oceania", "Siberia", 
#                                                 "Americas"), ordered = T)
panel$Pop <- factor(panel$Pop, levels = c("Anzick1", "Yaghan", "Aymara", "Quechua", "Wichi", "Toba",
                      "Yukpa", "Chane", "Guarani", "Ticuna", "Palikur",
                      "LagoaSanta", "Jamamadi", "Wayuu", 
                      "Guahibo",  "Bribri", "Zapotec1", "Tepehuano",
                      "Huichol", "Mixe", "Maya", "Aconcagua", "Ayayema",
                      "AncientYamana", "AncientKaweskar", "Taino", "USR1", 
                      "SpiritCave", "SouthWesternOntario", "PuntaSantaAna",
                      "TrailCreek", "Lovelock", "Botocudo","Piapoco","Pima", "Kogi", 
                      "Embera", "Waunana", "Teribe", "Guaymi",
                      "Maleku", "Cabecar","Karitiana", "Surui"), ordered = T)
# order botocudos by genome coverage
boto <- read.table("~/Projects/Botocudos/Files/Summaries/2019_02_25/Botocudos_summary_2019-02-25.table",
                   header = T)
boto$ID <- boto$Target
panel <- join(panel, boto[, c("ID", "hits_coverage_nuclear")], by = "ID")
panel$hits_coverage_nuclear[is.na(panel$hits_coverage_nuclear)] <- 100

pop_order <- order(panel$hits_coverage_nuclear, panel$Region, panel$Pop, decreasing = F)
names(pop_order) <- panel$Pop[pop_order] 
region_names <- tapply(1:nrow(panel), panel$Pop[pop_order], mean)
pop_names <- tapply(1:nrow(panel), panel$Pop[pop_order], mean)
pop_line <- tapply(1:nrow(panel), panel$Pop[pop_order], max)
panel <- panel[order(panel$Pop, panel$Pop),]

ten_colors <- c("#ED8F47", "#3F8EAA", 
                "#E52829","#79C360",
                "#9471B4","#FDB762", "#A6CEE3","#DDD399",
                "#B89B74", "#B15928")

Colors <- list("2_2" = c("#ED8F47", "#9471B4"), #orange,purple
               "3_14" = c("#ED8F47", "#79C360", "#9471B4"),
               #orange,purple,green
               "4_2" = c("#ED8F47", "#3F8EAA", "#9471B4","#79C360"),
               #orange,blue,green,purple
               "5_6" = c("#ED8F47",  "#3F8EAA",
                         "#9471B4", "#E52829","#79C360"),
               #orange,blue,red,green,purple
               "6_19" = c("#ED8F47", "#3F8EAA",
                         "#FDB762", "#E52829", "#79C360","#9471B4"),
               #orange,blue,purple,red,green,melon
               "7_6" =  c("#ED8F47", "#3F8EAA", "#79C360", 
                           "#E52829",  "#9471B4", "#FDB762",
                           "#A6CEE3"), 
               #orange,blue,green,red,purple,melon,lightblue
               "8_17" = c("#ED8F47",  "#3F8EAA", "#79C360",
                          "#E52829",
                          "#9471B4", "#FDB762", "#A6CEE3", "#DDD399"),
               #orange,blue,green,red,purple,melon,lightblue,beige
               "9_7" = c("#ED8F47",  "#3F8EAA", 
                         "#9471B4","#79C360","#E52829","#FDB762",
                         "#A6CEE3", "#DDD399", "#B89B74"),
               #orange,blue,purple,green,red,melon,lightblue,beige,farkbeige
               "10_3" = c("#ED8F47", "#3F8EAA", 
                           "#E52829","#79C360",
                           "#9471B4","#FDB762", "#A6CEE3","#DDD399",
                           "#B89B74", "#B15928")
               
)

png("~/Projects/Botocudos/Plots/ngsadmix/24boto_Maanasa_americas_fromVictor_americas.png",
    width = 27, height = 19, res = 300, units = "in")
par(mfrow = c(10,1))
for(k in paste(seq(2, 10 ), c( 27,85,83,40,94,77,86,14,60 ), sep = '_')){
  par(mar = c(3, .5, 3, 2))
  name <- paste(sufix, k, ".qopt", sep = "")
  
  
  admix<-as.matrix(read.table(name))#[pop_order,]
  #meanComponents <- colMeans(admix)
  
  admix <- admix[,order(meanComponents, decreasing = T)]
  #admix <- admix
  admix <- t(admix)
  barplot(admix[,1:nInd], space=0,border=NA, col = Colors[[k]], ylab=NULL,
          main = paste("k =", sub("_.*", "",k)), horiz = F, cex.main = 3,
          cex.lab = 2)
  abline(v = pop_line,
         col = "white")
}
#barplot(admix[,1:nInd], space=0,border=NA, col = c("Red", "Blue"))

text(region_names, -0.5, names(pop_names), xpd=NA, cex = 2, srt = 90)
dev.off()

text(pop_names, -0.5, names(pop_names), xpd=NA, cex = 2, srt = 90)

#-----------------------------------------------------------------------------#
library(ggplot2)
library(plyr)
library(adegenet)
#------------------------------------------------------------------------------#
# Read in files for panel
maanasa <- read.table("~/Projects/Botocudos/Files/ADMIX/2019_03_17/Worldwide/Subset_names.txt")
boto <- read.table("~/Projects/Botocudos/Files/ADMIX/2019_03_17/Worldwide/Botocudos_fromVictor.bam.list")
victor_panel <- read.table("~/Projects/Botocudos/Files/ADMIX/2019_03_17/ind_pop_region.txt")
maanasa_panel <- read.table("~/Projects/Botocudos/Files/ADMIX/2019_03_17/Maanasa_mask1_flip.panel", header =T)
moreno <- read.table("~/Projects/Botocudos/Files/ADMIX/2019_03_17/Maanasa.americas.Moreno_names.txt")

colnames(victor_panel) <- c("ID", "Pop", "Region")
colnames(moreno) <- c("Pop", "Language")
colnames(maanasa_panel) <- c("ID", "Fam", "Pop", "Region")
maanasa_panel$ID <- paste(maanasa_panel$ID, maanasa_panel$ID, sep = "_")
colnames(boto) <- "ID"
colnames(maanasa) <- "ID"
boto$ID <- sub(".*/", "", boto$ID)
boto$ID <- sub(".bam", "", boto$ID)

M <- join(maanasa, maanasa_panel, by = "ID")
M$Language <- M$Region
M$Fam <- NULL
B <- join(boto, victor_panel, by = "ID")
levels(B$Pop) <- c(levels(B$Pop), "Botocudo")
B[is.na(B$Pop), "Pop"] <- "Botocudo"
levels(B$Region) <- c(levels(B$Region), "Botocudo")
B[is.na(B$Region), "Region"] <- "Botocudo"
B <- join(B, moreno, by = "Pop")
levels(B$Language) <- c(levels(B$Language), levels(B$Pop))
B$Language[is.na(B$Language)] <- B$Pop[is.na(B$Language)]
panel <- rbind(B, M)
panel$Language <- factor(panel$Language,
                         levels = c("Botocudo", "LagoaSanta", "Taino",
                                    "EcuatorialTucanoan", "Andean",
                                    "CentralAmerind", "NorthernAmerind",
                                    "Americas", "Aconcagua", "Yukpa","AncientKaweskar",
                                    "AncientYamana", "Ayayema", "Anzick1", "Lovelock",
                                    "TrailCreek", "SpiritCave", "PuntaSantaAna",
                                    "SouthWesternOntario",
                                    "USR1", "Siberia", "Han","EastAsia", "SoutheastAsia",
                                    "SouthAsia", "Papuan", "Andaman", 
                                    "Australian", "Oceania","French",
                                    "Europe", "Yoruba", "Africa"), ordered = T)
panel$Pop <- factor(panel$Pop, levels = c("Botocudo", "LagoaSanta", "Karitiana",
                                          "Surui", "Taino","Chane", "Piapoco",
                                          "Aymara", "Quechua", "PuntaSantaAna", "Guahibo", "Palikur",
                                          "Guaymi", "Cabecar", "Teribe", "Bribri",
                                          "Maleku", "Ticuna", "Kogi", "Embera", "Waunana",
                                          "Wayuu", "Yaghan", "Jamamadi", "Guarani", "Wichi",
                                          "Toba", "Aconcagua", "Yukpa", "Ayayema",
                                          "Maya",
                                          "Zapotec1","Tepehuano", "Mixe",
                                          "Huichol", "Pima", "AncientKaweskar", 
                                          "AncientYamana", "Anzick1", "Lovelock", "TrailCreek", 
                                          "SpiritCave", "SouthWesternOntario", "USR1",
                                          "Koryaks", "Han", "Japanese", "Aeta", 
                                          "Pathan", "French","Papuan",
                                          "Papuans",  "Andaman", 
                                          "Australian", "Yoruba", "Yorubas"),
                    ordered = T)
panel$Region <- factor(panel$Region, 
                       levels = c("Botocudo", "Americas", "Siberia",
                                  "EastAsia", "SoutheastAsia", "SouthAsia",
                                  "Europe",
                                  "Oceania", "Africa"), ordered = T)
# Widths
panel$width <- 1.2
panel$width[panel$Region != "Americas"] <- 2
panel$width[panel$Pop == "Botocudo"] <- 4

# Remove ugly bars
toremove <- which(panel$ID == "Papuan" |panel$ID == "Han" |panel$ID == "Yoruba" |
                    panel$Region == "SouthAsia" | panel$Region == "SoutheastAsia" |
                    panel$Region == "Siberia" |panel$ID == "Australian")
panel <- panel[-toremove,]

#------------------------------------------------------------------------------#
# Order components
order_comp <- function(admix, regions){
  Ord <- c()
  for(r in regions){
    meanComp <- apply(as.matrix(admix[,1:ncomp]), 2,function(x) mean(x[admix$Region == r]))
    comp <- which(meanComp == max(meanComp))
    if(! comp %in% Ord){
      Ord <- c(Ord, comp)
    }
  } 
  return(Ord)
}
#------------------------------------------------------------------------------#
setwd("~/Projects/Botocudos/Files/ADMIX/2019_03_17/Worldwide/")
boto <- read.table("~/Projects/Botocudos/Files/Summaries/2019_02_25/Botocudos_summary_2019-02-25.table", header = T)
boto$ID <- boto$Target
boto <- boto[order(boto$hits_coverage_nuclear),]
sufix <- "Subset_Botocudos_fromVictor_k"
#---------------------------------------------------
admix <- read.table(name)
nInd <- dim(admix)[1]
ncomp <-  dim(admix)[2]
colnames(admix) <- seq(1,ncomp)
admix <- cbind(admix, panel)
admix <- admix[order(admix$Region, admix$Language, admix$Pop),]

barplot(t(as.matrix(admix[,1:ncomp])), space=0, border = NA)

# Colors
ten_colors <- c("#ED8F47", "#3F8EAA", 
                "#E52829","#79C360",
                "#9471B4","#FDB762", "#A6CEE3","#DDD399",
                "#B89B74", "#B15928")




ten_colors <- funky(5)
png("~/Desktop/admixture.png", width = 20, height = 20, res = 300, units = "in")
par(mfrow = c(5,1))
for(k in paste(seq(2,5), c( 27,85,83,40)#,94,77,86, 14,60,12,12,1,12,12,1,9,16,12,6 )
               , sep = '_')){
  par(mar = c(3, .5, 3, 2))
  name <- paste(sufix, k, ".qopt", sep = "")
  admix <- read.table(name)
  nInd <- dim(admix)[1]
  ncomp <-  dim(admix)[2]
  colnames(admix) <- seq(1,ncomp)
  admix <- admix[-toremove,]
  admix <- cbind(admix, panel)
  admix <- admix[order(admix$Pop, admix$Language, admix$Region, decreasing = F),]
  admix <- admix[-which(admix$ID %in% head(boto$ID, 7)),]
  #meanComponents <- apply(as.matrix(admix[,1:ncomp]), 2, function(x) median(x))#colMeans(as.matrix(admix[,1:ncomp]))
  admix[, 1:ncomp] <- admix[, order_comp(admix, c("Botocudo", "Africa", "Oceania", "EastAsia", "Europe"))]#admix[, order(meanComponents, decreasing = T)]
  admix <- join(admix, boto[, c("ID", "hits_coverage_nuclear")], by = "ID")
  admix$hits_coverage_nuclear[is.na(admix$hits_coverage_nuclear)] <- 100
  admix <- admix[ order(admix$hits_coverage_nuclear),]
  pop_line <- sapply(unique(admix$Region), 
                     function(X) which(
                       admix$ID == 
                         tail(admix$ID[admix$Region == X], 1)))
  names(pop_line) <- unique(admix$Region)
  pop_line <- pop_line[order(pop_line)]
  w <- sapply(names(pop_line), function(x) sum(admix$width[admix$Region == x]))
  pop_line <- cumsum(pop_line*w)
  barplot(as.matrix(t(admix[1:nInd,1:ncomp])), space=0,border=NA, 
          xlab=NULL,
          ylab = paste("k =", sub("_.*", "",k)), horiz = F, cex.main = 3,
          cex.lab = 2, col = ten_colors, axes = F, axisnames= F,
          width = admix$width)
  abline(v = cumsum(w),
         col = "white")
}

rownames(admix) <- 1:dim(admix)[1]
pop_names <- sapply(unique(admix$Region), 
                    function(X) mean(as.numeric(rownames(admix[admix$Region == X,]))))
names(pop_names) <- unique(admix$Region)
text(cumsum(w)-0.5*w, -0., names(pop_names), xpd=NA, cex = 1, srt = 90)
dev.off()
#text(pop_names, -0.5, names(pop_names), xpd=NA, cex = 1, srt = 90)

#------------------------------------------------------------------------------#
# K = 5
k <- "5_40"
panel$width <- 0.5
panel$width[panel$Region != "Americas"] <- 2
panel$width[panel$Pop == "Botocudo"] <- 4

 png("~/Projects/Botocudos/Plots/ngsadmix/17boto_Maanasa_k5.png",
     width = 6.5, height = 2.5, res = 300, units = "in", bg = NA)
par(mar = c(6, .5, 3, 2))
name <- paste(sufix, k, ".qopt", sep = "")
admix <- read.table(name)
nInd <- dim(admix)[1]
ncomp <-  dim(admix)[2]
colnames(admix) <- seq(1,ncomp)
admix <- admix[-toremove,]
admix <- cbind(admix, panel)
admix <- admix[order(admix$Pop, admix$Language, admix$Region, decreasing = F),]
admix <- admix[-which(admix$ID %in% head(boto$ID, 7)),]
meanComponents <- colMeans(as.matrix(admix[,1:ncomp]))
admix[, 1:ncomp] <- admix[, order(meanComponents, decreasing = T)]
admix[, 1:ncomp] <- admix[, order_comp(admix, (c("Botocudo", "EastAsia","Europe", "Oceania", "Africa")))]#admix[, order(meanComponents, decreasing = T)]

admix <- join(admix, boto[, c("ID", "hits_coverage_nuclear")], by = "ID")
admix$hits_coverage_nuclear[is.na(admix$hits_coverage_nuclear)] <- 100
admix <- admix[ order(admix$hits_coverage_nuclear),]
pop_line <- sapply(unique(admix$Region), 
                   function(X) which(
                     admix$ID == 
                       tail(admix$ID[admix$Region == X], 1)))
names(pop_line) <- unique(admix$Region)
pop_line <- pop_line[order(pop_line)]
w <- sapply(names(pop_line), function(x) sum(admix$width[admix$Region == x]))
pop_line <- cumsum(pop_line*w)
barplot(as.matrix(t(admix[1:nInd,1:ncomp])), space=0,border=NA, 
        xlab=NULL,
        ylab = paste("k =", sub("_.*", "",k)), horiz = F, cex.main = 3,
        cex.lab = 2, col = brewer_pal(palette = "Set2")(5), axes = F, axisnames= F,
        width = admix$width)
abline(v = cumsum(w),
       col = "white")
rownames(admix) <- 1:dim(admix)[1]
pop_names <- sapply(unique(admix$Region), 
                    function(X) mean(as.numeric(rownames(admix[admix$Region == X,]))))
names(pop_names) <- unique(admix$Region)
text(cumsum(w)-0.5*w, -0.3, names(pop_names), xpd=NA, cex = 1)
 dev.off()
