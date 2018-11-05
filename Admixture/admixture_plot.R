# Admixture plot
admixture_plot <- function(){
  setwd("~/Projects/Botocudos/Files/Quack/2018_10_03/")
  k <- "2_10"
  selected <- c("Quack", "MA2776", "MA2777")#, "Bot15", "Bot17")
  sufix <- "Jorde_Wollstein_hg19_reheaded_Quack_Bot15_Bot17_k"
  n <- 3
  N <- 583
  
  name <- paste(sufix, k, ".qopt", sep = "")
  admix<-t(as.matrix(read.table(name)))
  colors<- funky(gsub("_*", "",k))
  #panel_names <- "90ind_thesis.txt"
  panel <- read.table("~/Projects/Botocudos/Files/MDS/2017_12_06/Wollstein/panel_Jorde_Wollstein.txt",
                      header = T)
  
  pop <- panel[27:dim(panel)[1], c("region", "indID", "population")]
  pop <- rbind(data.frame(region = selected, 
                          indID = selected,
                          population = selected),
               pop)
  pop_order <- factor(pop[,3], levels=c("Bambaran", "Dogon", 
                                        "YRI",
                                        "CEU", 
                                        "Slovenian",
                                        "Iraqi",
                                        "Pakistanis", "Kyrgyzstan",
                                        "Mongolian", "JPT",
                                        "CHB", "BOR",
                                        "Tongan_Samoan", "POL", "FIJ",
                                        "Totonac", "Bolivian", "NGH", "Nepalese",
                                        "Thai",
                                        selected))
  pop <- pop[order(pop_order),]
  pop$population <- factor(pop$population, levels = levels(pop_order))
  
  
  
  colors<- funky(7)
  
  Colors <- list("2_90" = c("#ED8F47", "#9471B4"), #orange,purple
                 "3_1" = c("#ED8F47", "#9471B4", "#79C360"),
                 #orange,purple,green
                 "4_90" = c("#ED8F47", "#3F8EAA","#79C360", "#9471B4"),
                 #orange,blue,green,purple
                 "5_63" = c("#ED8F47",  "#3F8EAA","#E52829","#79C360",
                            "#9471B4"),
                 #orange,blue,red,green,purple
                 "6_1" = c("#ED8F47", "#3F8EAA",
                           "#FDB762", "#E52829", "#79C360","#9471B4"),
                 #orange,blue,purple,red,green,melon
                 "7_96" =  c("#ED8F47", "#3F8EAA", "#79C360", 
                             "#E52829",  "#9471B4", "#FDB762",
                             "#A6CEE3"), 
                 #orange,blue,green,red,purple,melon,lightblue
                 "8_14" = c("#ED8F47",  "#3F8EAA", "#79C360",
                            "#E52829",
                            "#9471B4", "#FDB762", "#A6CEE3", "#DDD399"),
                 #orange,blue,green,red,purple,melon,lightblue,beige
                 "9_2" = c("#ED8F47",  "#3F8EAA", 
                           "#9471B4","#79C360","#E52829","#FDB762",
                           "#A6CEE3", "#DDD399", "#B89B74"),
                 #orange,blue,purple,green,red,melon,lightblue,beige,farkbeige
                 "10_21" = c("#ED8F47", "#3F8EAA", 
                             "#E52829","#79C360",
                             "#9471B4","#FDB762", "#A6CEE3","#DDD399",
                             "#B89B74", "#B15928")
                 
  )
  
  
  
  par(mfrow = c(2,1))
  
  for(k in paste(7,77, sep = '_')){
    par(mar = c(7, 0.5, 3, 2))
    name <- paste(sufix, k, ".qopt", sep = "")
    
    admix<-as.matrix(read.table(name))
    admix <- admix[order(pop_order),]
    admix <- t(admix)
    indexes <- col2pop(colors, pop, N+n, 10)
    barplot(admix[indexes[1:dim(admix)[1]],],
            width = c(rep(.1,N),rep(3,n)),
            col=funky(7),space=0,border=NA,
            ylab=NULL,
            main = paste("k =", sub("_.*", "",k)), horiz = F, cex.main = 3,
            cex.lab = 2)
    
    pop_names <- 0.1*tapply(1:nrow(pop),pop[,3],mean)
    
    abline(v = 0.1*tapply(2:(nrow(pop)-n),pop[2:(nrow(pop)-n),3],max),
           col = "white")
    abline(v = 3*seq(1,n)+N*0.1, col = "white")
    
    
  }
  pop_names[match(tail(pop_names, n), pop_names)] <- tail(pop_names, n) + 2.5*(seq(1,n))
  text(pop_names, -0.5,names(pop_names),xpd=NA,cex = 2, srt = 90)
}