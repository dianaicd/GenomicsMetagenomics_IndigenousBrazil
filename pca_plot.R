pca_plot <- function(ind, plink=F, d, pcx = 1, pcy = 2, colors = F, label = F, panel_name="merged_noduplicates_reheaded"){
  
  panel <- read.table("~/Projects/Botocudos/Files/PCA/2017_09_15/panel_botos.txt",
                      header = T)
  
  popinfo <- read.csv("~/Projects/Botocudos/Files/PCA/2017_09_15/panel_botos.txt",
                      sep="\t", header=T)
  path <- paste(d, ind, sep = "")
  
  if(plink){
    eval <- read.table(paste(path, ".eigenval", sep = ""))
    evec <- read.table(paste(path, ".eigenvec", sep = ""))
    # Modifications due to plink not acceptind underscores
    # in the vcf with new header

    evec$V1 <- sub("-", "_", evec$V1)
    evec$V1 <- paste(evec$V1, evec$V2, sep = "_")
    evec$V2 <- NULL
  }else{
    eval <- read.table(paste(path, ".hg19_nuc.realigned.filtSites_",
                             panel_name,"_MERGED.eval", sep = ""))
    evec <- read.table(paste(path, ".hg19_nuc.realigned.filtSites_",
                             panel_name,"_MERGED.evec", sep = ""))
    # Modifications due to plink not acceptind underscores
    # in the vcf with new header
    evec$V1 <- sub("S-", "S_", evec$V1)
    evec$V1 <- sub("-", "_", evec$V1)
    evec$V1 <- gsub(":", "_", evec$V1)
    evec$V1 <- sub("_1", "-1", evec$V1)
    evec$V1 <- sub("_2", "-2", evec$V1)
    evec$V1 <- sub("_3", "-3", evec$V1)
  }
  
  n <- dim(evec)[2]
  pca <- evec[,seq(2,n)]
  
  individual.name <- gsub("_MA.*", "", evec[,1])
  individual.name <- gsub(".*_HGDP", "HGDP", individual.name)
  individual.name <- gsub(".*_HG", "HG", individual.name)
  individual.name <- gsub(".*_NA", "NA", individual.name)
  individual.name <- gsub(".*_S_", "S_", individual.name)
  individual.name <- gsub(".*_B_", "B_", individual.name)
  individual.name <- gsub("0_HG", "HG", individual.name)
  
  evec[,1] <- individual.name
  sample.idx <- match(ind, gsub(".*_","",evec[,1]))
  
  
  colnames(evec)<- c("indID", paste("pca", seq(1,n-1), sep=""))

  evec <- join(evec, popinfo, by = "indID")
  regions <- as.character(unique(evec$region))
  populations <- as.character(unique(evec$population))
  
  botox <- evec[sample.idx, pcx+1]
  botoy <- evec[sample.idx, pcy+1]
  
  if(!colors){
    colors <- funky(length(regions)-1)
    colors <- c("black", colors)
  }

  if(label){
    p <- ggplot(evec[rev(rownames(evec)),], aes(x = evec[rev(rownames(evec)),pcx+1], 
                                                y = evec[rev(rownames(evec)),pcy+1], 
                          color = region, alpha = region, label = population)) +
      geom_text() +
      scale_color_manual(values = colors, name = "Region",
                         breaks = regions,
                         labels = regions) +
      scale_alpha_manual(name = "Region",
                         values = c(1, rep(0.7, length(evec$region) - 1)),
                         breaks = regions,
                         labels = regions) +
      labs(title = ind, x = pct_var_pc(eval, pcx),
           y = pct_var_pc(eval, pcy))
  }else{
  
  # evec$region[is.na(evec$region)] <- "Botocudos"
  # index <- grepl("Boto", unique(evec$region))
  # evec$region <- factor(evec$region, levels = c("Botocudos", 
  #                                               as.character(unique(evec$region)[!index]))[1:11])
  p <- ggplot(evec, aes(x = evec[,pcx+1], 
                                              y = evec[,pcy+1], 
                        color = region, alpha = region)) +
    geom_point(size = 4) +
    scale_color_manual(values = colors, name = "Region",
                       breaks = regions,
                       labels = regions) +
    scale_alpha_manual(name = "Region",
                       values = c(1, rep(0.7, length(evec$region) - 1)),
                                  breaks = regions,
                                  labels = regions) +
    labs(title = ind, x = pct_var_pc(eval, pcx),
          y = pct_var_pc(eval, pcy))  + 
    geom_point(aes(x = botox, y = botoy), color = "black", size=4, show.legend = F)#+
    # geom_point(aes(x = botox, y = botoy, color = "black"))
  }
  return(p)
  
}
