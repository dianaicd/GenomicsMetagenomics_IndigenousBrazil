library(data.table)
library(gtools)
library(FactoMineR)

setwd("~/Vital-IT/sapfo/australia/")

## read file
data<-as.data.frame(t(read.table("flagstat.txt.bk", sep="\t", row.names=1, header=T)))
data<-cbind(data, "group"=gsub('[0-9|a-z|.]+', '', rownames(data)))
data<-data[mixedorder(rownames(data)),]

group<-data.table("group"=data$group, "idx"=1:nrow(data))
groups<-group[,list(idx=mean(idx)), by='group']

#data<-as.data.table(data)

pdf("Austrlaia_mappings.bk.pdf", paper="a4", width=0, height=0)
par(mfrow = c(3,1)) #add room for the rotated labels

m<-barplot(data[,1], col=data$group, main="total number of reads")
axis(1, at=m[groups$idx], labels=groups$group, las=3)

barplot(1-data[,5]/data[,1], col=data$group, main="% of unmapped reads")
axis(1, at=m[groups$idx], labels=groups$group, las=3)

barplot(data[,1]-data[,5], col=as.factor(data[,14]), main="number of unmapped reads")
axis(1, at=m[groups$idx], labels=groups$group, las=3)

dev.off()

range(1-data[,5]/data[,1])
mean(data[,5]/data[,1])

range(data[,1]-data[,5])
mean(data[,1]-data[,5])

m<-barplot(data[,5], col=data$group, main="number of mapped reads")
axis(1, at=m[groups$idx], labels=groups$group, las=3)


pdf("Australia_mappings_flagstat.pdf", paper="a4", width=0, height=0)
par(mfrow = c(3,1)) #add room for the rotated labels
for(i in 1:(ncol(data)-1)){
  m<-barplot(data[,i], col=data$group, main=colnames(data[i]))
  axis(1, at=m[groups$idx], labels=groups$group, las=3)
}
dev.off()


############################################################################
## get number of classifeid and unclassfied reads
unclassified<-vector(length=length(all.list))
classified<-vector(length=length(all.list))
for(i in 1:length(all.list)){unclassified[i]<-all.list[[i]]$reads[1]; classified[i]<-all.list[[i]]$reads[2]}


## average number of unmapped reads
head(classified+unclassified,1)

## average number of unclssified reads
head(unclassified,1)

## average total number of reads
head(data$`in total `,1)

## average total number of mapped reads
(head(data$`in total `,1) - head(data$mapped,1))/2



############################################################################
############################################################################
## read centrifuge summary output
files<-list.files(path="~/Vital-IT/sapfo/australia/centrifuge_summary", pattern="*.summary", full.names=T, recursive=FALSE)
inds<-gsub(".summary","",gsub("/Users/sneuensc/Vital-IT/sapfo/australia/centrifuge_summary/centrifuge_", "", files))

col=6
m<-NULL
for(i in 1:length(files)){
  print(inds[i])
  #cen<-read.table(files[i], h=T, sep="\t")
  cen<-fread(files[i],sep="\t")[,c(1:4,col), with=F]
  setnames(cen, c(colnames(cen)[1:4],inds[i]))
  #cen<-fread(files[i])
  #if(i==1) {m<-cen[,c(1:5)]} else {m<-merge(m, cen[,c(2,5)], all=T)}
  if(i==1) {m<-cen} else {m<-merge(m, cen, all=T)}
  #taxRank<-cen[,list(.N), by=taxRank]}
#  m<-cbind(m,taxRank$N)
}
m[is.na(m)]<-0



## get the number of each taxRank
pdf("Australia_taxRank.pdf", paper="a4", width=0, height=0)
par(mfrow = c(3,1)) #add room for the rotated labels

(taxRank<-m[,-c(1,2,4), with=F][,lapply(.SD, function(x)return(sum(x>0))), by=taxRank])

for(i in 1:nrow(taxRank)){
  m<-barplot(as.numeric(taxRank[i,-1, with=F]), col=as.factor(data[,14]), main=taxRank[i,1])
  axis(1, at=m[groups$idx], labels=groups$group, las=3)
}
dev.off()
  

############################################################################
############################################################################
## read centrifuge kraken report

files<-list.files(path="~/Vital-IT/sapfo/australia/centrifuge_report", pattern="*.report", full.names=T, recursive=FALSE)
inds<-gsub(".report","",gsub("/Users/sneuensc/Vital-IT/sapfo/australia/centrifuge_report/centrifuge_", "", files))

taxRanks<-matrix(c("U", "unclassified", 
                "D", "domain",
                "P", "pyhlum",
                "C", "class",
                "F", "family",
                "G", "genus",
                "S", "species"), byrow=T, ncol=2)
absolute<-2


pdf("Australia_taxRank.pdf", paper="a4", width=0, height=0)
par(mfrow = c(3,1)) #add room for the rotated labels

for(r in 1:nrow(taxRanks)){
  for(i in 1:length(files)){
    print(inds[i])
    #cen<-read.table(files[i], h=T, sep="\t")
    cen<-fread(files[i],sep="\t")
    setnames(cen, c("prop","number","V1","taxRank","taxID","taxon"))
    cen<-cen[taxRank==taxRanks[r,1]]
    cen<-cen[,c(5,6,absolute), with=F]
    setnames(cen, c("taxID","taxon", inds[i]))
    if(i==1) {m<-cen} else {m<-merge(m, cen, all=T)}
    #taxRank<-cen[,list(.N), by=taxRank]}
    #  m<-cbind(m,taxRank$N)
  }
  my.row<-m$taxon
  m<-m[,-(1:2), with=F]
  m<-m[,mixedorder(colnames(m)), with=F]  
  m[is.na(m)]<-0
  group<-data.table("group"=gsub('[0-9|a-z|.|-]+', '', colnames(m)), "idx"=1:ncol(m))
  groups<-group[,list(idx=mean(idx), .N), by='group']

  #(taxRank<-m[,-c(1,2,4), with=F][,lapply(.SD, function(x)return(sum(x>0))), by=taxRank])
  a<-unlist(lapply(m, function(x)return(sum(x>0))))
  q<-barplot(as.vector(a), col=as.factor(data[,14]), main=taxRanks[r,2])
  axis(1, at=q[groups$idx], labels=groups$group, las=3)
}
dev.off()

colours <- c("red", "orange", "blue", "yellow", "green")
q<-barplot(prop.table(as.matrix(m),2), col=palette(), axisnames=F, main=taxRanks[r,2], ylab="proportion of reads")
bin<-(q[2]-q[1])
groups$Ncum<-cumsum(groups$N)
f1<-q[1]-bin/2
bins<-c(f1, groups$Ncum*bin+f1)
axis(1, at=bins, labels=rep("",length(bins)))
axis(1, at=q[groups$idx], labels=groups$group, las=3, tick=F)
abline(v=bins, col="black", lwd=2)
legend(0,1 , my.row, cex=1.3, bty="n", fill=palette())

############################################################################
## get the taxonomy tree for each entry
files<-list.files(path="~/Vital-IT/sapfo/australia/centrifuge_report", pattern="*.report", full.names=T, recursive=FALSE)
inds<-gsub(".report","",gsub("/Users/sneuensc/Vital-IT/sapfo/australia/centrifuge_report/centrifuge_", "", files))
taxRanks<-matrix(c("U", "unclassified",
                   "D", "domain",
                   "K", "kingdom",
                   "P", "phylum",
                   "C", "class",
                   "O", "order",
                   "F", "family",
                   "G", "genus",
                   "S", "species"), byrow=T, ncol=2)

# all.list<-NULL
# for(f in 1:length(files)){
#   print(inds[f])
#   cen<-fread(files[f],sep="\t")
#   setnames(cen, c("prop","number","new","taxRank","taxID","taxon"))
#   cen<-cen[taxRank!="U" & taxRank!="-",]
#   m<-NULL
#   taxas<-NULL
#   for(i in 1:nrow(cen)){
#     tax=cen[i]$taxRank
#     rank<-which(taxRanks[,1]==tax)
#     
#     ## adjust taxas vector
#     while(length(taxas)>0 && 
#           which(unlist(strsplit(taxas[length(taxas)],split="_"))[1]==taxRanks[,1]) >= rank){
#       taxas<-taxas[-length(taxas)]
#     }
#     taxas<-c(taxas, paste(tax,cen[i]$taxon, sep="_"))
#     
#     cen$taxonomy[i]<-paste(taxas, collapse="; ")
#   }
#   all.list[[inds[f]]]<-cen
# }
# saveRDS(all.list, "all_centrifuge.rds")


all.list2<-NULL
for(f in 1:length(files)){
  print(inds[f])
  all.list2[[inds[f]]]<-read_report2(files[f], collapse = TRUE, keep_ranks = c("D", "K", "P", "C", "O", "F", "G", "S"), 
               min.depth = 0, filter_taxon = NULL,
               has_header = NULL, add_rank_columns = T)
}
saveRDS(all.list2, "all_centrifuge2.rds")

m[m]

############################################################################

all.list<-readRDS("all_centrifuge2.rds")
all.list2<-readRDS("all_centrifuge2.rds")

pdf("Australia_taxRank2.pdf", paper="a4", width=0, height=0)
par(mfrow = c(3,1)) #add room for the rotated labels

dev.off()


# get.centrifuge<-function(all.list, filter.taxRank=NULL, filter.taxon=NULL, filter.taxonomy=NULL, 
#                          invert.taxRank=F, invert.taxon=F, invert.taxonomy=F, absolute=2)
# {
#   m<-NULL
#   for(r in 1:nrow(taxRanks)){
#     r<-2
#     for(i in 1:length(all.list)){
#       print(inds[i])
#       cen<-all.list[[inds[i]]]
#       if(!is.null(filter.taxRank)) cen<-cen[grep(paste0(filter.taxRank, collapse="|"), cen$taxRank, invert=invert.taxRank),]
#       if(!is.null(filter.taxon)) cen<-cen[grep(paste0(filter.taxon, collapse="|"), cen$taxon, invert=invert.taxon),]
#       if(!is.null(filter.taxonomy)) cen<-cen[grep(paste0(filter.taxonomy, collapse="|"), cen$taxonomy, invert=invert.taxonomy),]
#       cen<-cen[,c(5,6,absolute), with=F]
#       cen<-as.data.frame(cen)
#       cen$taxID<-unlist(cen$taxID)
#       cen<-cen[order(cen$taxID),]
#       setnames(cen, c("taxID","taxon", inds[i]))
#       if(i==1) {m<-cen} else {m<-merge(m, cen, all=T)}
#     }
#     my.leg<-m$taxon
#     m<-m[,-(1:2)]
#     m<-m[,mixedorder(colnames(m))]  
#     m[is.na(m)]<-0
#     m<-apply(m, 2,as.numeric)
#     group<-data.table("group"=gsub('[0-9|a-z|.|-]+', '', colnames(m)), "idx"=1:ncol(m))
#     groups<-group[,list(idx=mean(idx), .N), by='group']
#   }
#   return (list("my.leg"=my.leg, "m"=as.data.table(m), "groups"=groups))
# }

get.centrifuge<-function(all.list, filter.taxRank=NULL, filter.taxon=NULL, filter.taxonomy=NULL, 
                         remove.taxRank=NULL, remove.taxon=NULL, remove.taxonomy=NULL, absolute=2)
{
  m<-NULL
  for(i in 1:length(all.list)){
    print(inds[i])
    cen<-all.list2[[inds[i]]]
    
    ## get the taxa to remove
    rem<-NULL
    if(!is.null(remove.taxRank)) rem<-rbind(rem, cen[grep(paste0(remove.taxRank, collapse="|"), cen$rank, ignore.case=T),])
    if(!is.null(remove.taxon)) rem<-rbind(rem, cen[grep(paste0(remove.taxon, collapse="|"), cen$name, ignore.case=T),])
    if(!is.null(remove.taxonomy)) rem<-rbind(rem, cen[grep(paste0(remove.taxonomy, collapse="|"), cen$taxonstring, ignore.case=T),])
    
    ## taxa to fitler for
    if(!is.null(filter.taxRank)) cen<-cen[grep(paste0(filter.taxRank, collapse="|"), cen$rank, ignore.case=T),]
    if(!is.null(filter.taxon)) cen<-cen[grep(paste0(filter.taxon, collapse="|"), cen$name, ignore.case=T),]
    if(!is.null(filter.taxonomy)) cen<-cen[grep(paste0(filter.taxonomy, collapse="|"), cen$taxonstring, ignore.case=T),]
    
    ## remove the reads of the taxa to filter out
    if(!is.null(rem)){
      for(t in 1:nrow(rem)){ ## for each entry to remove
        l<-unlist(strsplit(rem$taxonstring[t], "|", fixed=T))
        nb.reads<-rem$reads_stay[t]
        for(ll in 1:length(l)){    # for each domain above it remove the reads
          idx<-which(cen$rank==toupper(substring(l[ll], 1, 1)) & substring(cen$name,3)==substring(l[ll], 3))
          if(length(idx)>0){
              cen$reads[idx]<-cen$reads[idx]-nb.reads
          }
        }
      }
    }
    cen <- cen[!(cen$reads<=0),]
    cen.example<-cen
    
    if(nrow(cen)==0) return (NULL)
    cen<-cen[,c(5,6,absolute)]
    cen<-cen[order(cen$taxonid),]
    cen$name<-unlist(lapply(cen$name, function(x) return (substring(x, 3))))
    setnames(cen, c("taxID","taxon", inds[i]))
    if(i==1) {m<-cen} else {m<-merge(m, cen, all=T)}
  }
  my.leg<-m$taxon
  m<-m[,-(1:2)]
  m<-m[,mixedorder(colnames(m))]  
  m[is.na(m)]<-0
  group<-data.table("group"=gsub('[0-9|a-z|.|-]+', '', colnames(m)), "idx"=1:ncol(m))
  groups<-group[,list(idx=mean(idx), .N), by='group']
  
  return (list("my.leg"=my.leg, "m"=as.data.table(m), "groups"=groups, "example"=cen.example))
}

get.centrifuge2<-function(all.list, filter.taxRank=NULL, filter.taxon=NULL, filter.taxonomy=NULL, 
                          remove.taxRank=NULL, remove.taxon=NULL, remove.taxonomy=NULL, absolute=2)
{
  q<-get.centrifuge(all.list, filter.taxRank, filter.taxon, filter.taxonomy, 
                    remove.taxRank, remove.taxon, remove.taxonomy, absolute)
  if(is.null(q)) return (NULL)
  cum.sum<-apply(q$m,1, sum)
  my.sort<-order(cum.sum, decreasing=T)
  my.legend<-paste0(q$my.leg, " ", format(100*cum.sum/sum(cum.sum, na.rm=T), digits=0, scientific=F),"%")[my.sort]
  m<-prop.table(as.matrix(q$m[my.sort,]),2)
  
  return (list("my.legend"=my.legend, "m"=m, "groups"=q$groups, "example"=q$example, "read.mean"=mean(cum.sum), "read.min"=min(cum.sum), "read.max"=max(cum.sum)))
}


q<-get.centrifuge(all.list2, filter.taxRank="D", remove.taxon=c("Homo sapiens", "synthetic construct", "phiX174", "Phix174"), absolute=2)
apply(q$m,1,mean)

###############################################
## get classified
png(paste0("classified_barplot.png"), width = 1024, height = 1024, res=100)
qw<-get.centrifuge(all.list2, filter.taxRank="U", absolute=1)
qq<-barplot(as.numeric(100-qw$m)/100, col=rep(palette(),2)[group$group], axisnames=F, main="Classified", ylab="proportion of reads", ylim=0:1)
bin<-(qq[2]-qq[1])
q$groups$Ncum<-cumsum(q$groups$N)
f1<-qq[1]-bin/2
bins<-c(f1, q$groups$Ncum*bin+f1)
axis(1, at=bins, labels=rep("",length(bins)))
axis(1, at=qq[q$groups$idx], labels=q$groups$group, las=3, tick=F)
abline(v=bins, col="black", lwd=2)
lines(rep(0,length(bins)))
dev.off()

###############################################
## get classified II: reads to filter out
png(paste0("reads_to_filter_out_barplot.png"), width = 1024, height = 1024, res=100)
unclass.abs<-get.centrifuge(all.list2, filter.taxRank="U", absolute=2)
unclass.rel<-get.centrifuge(all.list2, filter.taxRank="U", absolute=1)
tot.reads<-unclass.abs$m / unclass.rel$m * 100

filterout<-get.centrifuge(all.list2, filter.taxRank="S", filter.taxon=c("Homo sapiens", "synthetic construct", "phiX174", "Phix174"), absolute=2)

tot.usable<-tot.reads-unclass.abs$m-apply(filterout$m,2,sum)
qq<-barplot(as.numeric(tot.usable/tot.reads), col=rep(palette(),2)[group$group], axisnames=F, main="Classified & usable", ylab="proportion of reads", ylim=0:1)

qq<-barplot(as.numeric(tot.usable), col=rep(palette(),2)[group$group], axisnames=F, main="Classified & usable", ylab="total number of reads")
mean(as.numeric(tot.usable))
range(as.numeric(tot.usable))
      
for(i in 1:3){
  print(mean(as.numeric(filterout$m[i,]/tot.reads*100), na.rm=F))
}

apply(filterout$m,1,function(x)return(mean(x/tot.reads, na.rm=T))))
filterout$my.leg

                       
ssum/unclass.tot*100
qq<-barplot(as.numeric(rel.usable), col=rep(palette(),2)[group$group], axisnames=F, main="Reads to filter out", ylab="proportion of reads")
            qq<-barplot(qqqq, col=palette()[-1], axisnames=F, main="Reads to filter out", ylab="proportion of reads")
                        
qq<-barplot(q$m, col=palette()[-1], axisnames=F, main="Reads to filter out", ylab="proportion of reads")
bin<-(qq[2]-qq[1])
q$groups$Ncum<-cumsum(q$groups$N)
f1<-qq[1]-bin/2
bins<-c(f1, q$groups$Ncum*bin+f1)
axis(1, at=bins, labels=rep("",length(bins)))
axis(1, at=qq[q$groups$idx], labels=q$groups$group, las=3, tick=F)
abline(v=bins, col="black", lwd=2)
lines(rep(0,length(bins)))
par(oldpar)
dev.off()

png(paste0("reads_to_filter_out_legend.png"), width = 1024, height = 1024, res=200)
plot(1, type="n", axes=F, xlab="", ylab="")
legend("topleft", q$my.legend, cex=1, bty="n", fill=palette()[-1], title=paste0(format(q$read.mean, scientific=T, digits=2), " reads per sample on average:"), title.adj=0)
dev.off()

# library(pvclust)
# fit <- pvclust(as.data.frame(q$m), parallel=T)
# #pdf("hierarchical_clustering.pdf", paper="a4r", width=0, height=0)
# plot(fit) # dendogram with p values
# #dev.off()



###############################################
## plot domain
png(paste0("domain_distribution_barplot.png"), width = 1024, height = 1024, res=100)
q<-get.centrifuge2(all.list2, filter.taxRank="D", remove.taxon=c("Homo sapiens", "synthetic construct", "phiX174", "Phix174"), absolute=2)
qq<-barplot(q$m, col=palette()[-1], axisnames=F, main=taxRanks[2,2], ylab="proportion of reads")
bin<-(qq[2]-qq[1])
q$groups$Ncum<-cumsum(q$groups$N)
f1<-qq[1]-bin/2
bins<-c(f1, q$groups$Ncum*bin+f1)
axis(1, at=bins, labels=rep("",length(bins)))
axis(1, at=qq[q$groups$idx], labels=q$groups$group, las=3, tick=F)
abline(v=bins, col="black", lwd=2)
lines(rep(0,length(bins)))
dev.off()

png(paste0("domain_distribution_legend.png"), width = 1024, height = 1024, res=200)
plot(1, type="n", axes=F, xlab="", ylab="")
legend("topleft", q$my.legend, cex=1, bty="n", fill=palette()[-1], 
       title=paste0(format(q$read.mean, scientific=T, digits=2), " reads per sample on average:"), justify="left", title.adj=0)
dev.off()




## hierarchical clustering
# fit <- pvclust(as.data.frame(q$m), parallel=T)
# plot(fit) # dendogram with p values

## pca
# pca <- prcomp(t(q$m), center = TRUE, scale. = TRUE) 
# plot(pca, type="l")
# my.max<-max(abs(c(pca$x[, 1], pca$x[, 2])))
# oldpar<-par(pty="s")
# plot(pca$rotation[, 1], pca$rotation[, 2], col=as.factor(group$group), main = "PCA", xlab = "PC1", ylab = "PC2")
# plot(pca$x[, 1], pca$x[, 2], col=as.factor(group$group), main = "PCA", xlab = "PC1", ylab = "PC2")
# par(oldpar)
# biplot(pca, var.axes=F, col=c("black","white"))

# library(FactoMineR)
# oldpar<-par(pty="s"); plot(PCA(t(q$m), graph = F),col.ind=as.factor(group$group)); par(oldpar)


#############################################
## pdfs
domains<-c("Bacteria","Viruses","Archaea","Eukaryota", "Viroids" )
for(domain in domains){
  pdf(paste0(domain,"_distribution.pdf"), paper="a4", width=0, height=0)
  layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
  for(i in 3:nrow(taxRanks)){
    cur.rank<-taxRanks[i,1]   ## "U" "D" "K" "P" "C" "O" "F" "G" "S"
    q<-get.centrifuge2(all.list2, filter.taxRank=cur.rank, remove.taxon=c("Homo sapiens", "synthetic construct", "phiX174", "Phix174"), 
                       filter.taxonomy=domain, absolute=2)
    if(length(q$my.legend)>24) q$my.legend<-q$my.legend[1:24]
    
    if(is.null(q)) {
      print(paste("no data for", taxRanks[i,2]))
      next
    }
    
    qq<-barplot(q$m, col=palette()[-1], axisnames=F, main=taxRanks[taxRanks[,1]==cur.rank,2], ylab="proportion of reads", ylim=0:1)
    bin<-(qq[2]-qq[1])
    q$groups$Ncum<-cumsum(q$groups$N)
    f1<-qq[1]-bin/2
    bins<-c(f1, q$groups$Ncum*bin+f1)
    axis(1, at=bins, labels=rep("",length(bins)))
    axis(1, at=qq[q$groups$idx], labels=q$groups$group, las=3, tick=F)
    abline(v=bins, col="black", lwd=2)
    lines(rep(0,length(bins)))
    
    plot(1, type="n", axes=F, xlab="", ylab="")
    legend("topleft", q$my.legend[1:24], cex=1, bty="n", fill=palette()[-1], title=paste0(format(q$read.mean, scientific=T, digits=2), " reads per sample on average:"), title.adj=0)
    
    oldpar2<-par(pty="s"); 
    plot(PCA(t(q$m), graph = F),col.ind=as.factor(group$group), title=paste("PCA of",taxRanks[taxRanks[,1]==cur.rank,2])); 
    par(oldpar2)
  }
  dev.off()
}

#############################################
## pdfs 2
domains<-c("Bacteria","Viruses","Archaea","Eukaryota", "Viroids" )
for(domain in domains){
  pdf(paste0(domain,"_distribution2.pdf"), paper="a4", width=0, height=0)
  layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
  for(i in 3:nrow(taxRanks)){
    cur.rank<-taxRanks[i,1]   ## "U" "D" "K" "P" "C" "O" "F" "G" "S"
    q<-get.centrifuge2(all.list2, filter.taxRank=cur.rank, remove.taxon=c("Homo sapiens", "synthetic construct", "phiX174", "Phix174"), 
                       filter.taxonomy=domain, absolute=2)
    if(length(q$my.legend)>24) q$my.legend<-q$my.legend[1:24]
    
    if(is.null(q)) {
      print(paste("no data for", taxRanks[i,2]))
      next
    }
    
    nb.col=2
    q$m<-q$m[1:nb.col,]
    q$m<-q$m[,colnames(q$m)!="NA04932" & colnames(q$m)!="P2077"]
    q$groups<-q$groups[q$groups$group!="NA" & q$groups$group!="P",]

    q$my.legend<-q$my.legend[1:nb.col]
    q$my.legend[nb.col+1]<-paste0("others ", 100-sum(as.integer(gsub("%","",lapply(strsplit(q$my.legend[1:nb.col]," "), tail, n = 1L)))),"%")

    qq<-barplot(q$m, col=palette()[-1], axisnames=F, main=taxRanks[taxRanks[,1]==cur.rank,2], ylab="proportion of reads", ylim=0:1)
    bin<-(qq[2]-qq[1])
    q$groups$Ncum<-cumsum(q$groups$N)
    f1<-qq[1]-bin/2
    bins<-c(f1, q$groups$Ncum*bin+f1)
    axis(1, at=bins, labels=rep("",length(bins)))
    axis(1, at=qq[q$groups$idx], labels=q$groups$group, las=3, tick=F)
    abline(v=bins, col="black", lwd=2)
    lines(rep(0,length(bins)))
    
    plot(1, type="n", axes=F, xlab="", ylab="")
    col<-palette()[-1][1:nb.col]
    col[nb.col+1]<-"white"
    legend("topleft", q$my.legend, cex=1, bty="n", fill=col, title=paste0(format(q$read.mean, scientific=T, digits=2), " reads per sample on average:"), title.adj=0)
    
    oldpar2<-par(pty="s"); 
   # m2<-q$m[1:2,]
   # m2[2,]<-1-m2[1,]
  #  m2<-m2[,colnames(m2)!="NA04932" & colnames(m2)!="P2077"]
    # plot(PCA(t(m2), graph = F),col.ind=as.factor(group$group), title=paste("PCA of",taxRanks[taxRanks[,1]==cur.rank,2])); 
    plot(PCA(t(q$m), graph = F),col.ind=as.factor(group$group), title=paste("PCA of",taxRanks[taxRanks[,1]==cur.rank,2])); 
    par(oldpar2)
  }
  dev.off()
}


domains<-c("Bacteria","Viruses","Archaea","Eukaryota", "Viroids" )
for(domain in domains){
  dir.create(domain, showWarnings=F)
  for(i in 3:nrow(taxRanks)){
    cur.rank<-taxRanks[i,1]   ## "U" "D" "K" "P" "C" "O" "F" "G" "S"
    q<-get.centrifuge2(all.list2, filter.taxRank=cur.rank, remove.taxon=c("Homo sapiens", "synthetic construct", "phiX174", "Phix174"), 
                       filter.taxonomy=domain, absolute=2)
    if(length(q$my.legend)>24) q$my.legend<-q$my.legend[1:24]
    
    if(is.null(q)) {
      print(paste("no data for", taxRanks[i,2]))
      next
    }
    
    png(paste0(domain,"/",domain,"_distribution_",taxRanks[i,2],"_barplot.png"), width = 1024, height = 1024, res=100)
    qq<-barplot(q$m, col=palette()[-1], axisnames=F, main=taxRanks[taxRanks[,1]==cur.rank,2], ylab="proportion of reads", ylim=0:1)
    bin<-(qq[2]-qq[1])
    q$groups$Ncum<-cumsum(q$groups$N)
    f1<-qq[1]-bin/2
    bins<-c(f1, q$groups$Ncum*bin+f1)
    axis(1, at=bins, labels=rep("",length(bins)))
    axis(1, at=qq[q$groups$idx], labels=q$groups$group, las=3, tick=F)
    abline(v=bins, col="black", lwd=2)
    lines(rep(0,length(bins)))
    dev.off()
    
    png(paste0(domain,"/",domain,"_distribution_",taxRanks[i,2],"_legend.png"), width = 1024, height = 1024, res=200)
    plot(1, type="n", axes=F, xlab="", ylab="")
    legend("topleft", q$my.legend, cex=1, bty="n", fill=palette()[-1], title=paste0(format(q$read.mean, scientific=T, digits=2), " reads per sample on average:"), title.adj=0)
    dev.off()
    
    png(paste0(domain,"/",domain,"_distribution_",taxRanks[i,2],"_pca.png"), width = 1024, height = 1024, res=200)
    oldpar2<-par(pty="s"); 
    plot(PCA(t(q$m), graph = F),col.ind=as.factor(group$group), title=paste("PCA of",taxRanks[taxRanks[,1]==cur.rank,2])); 
    par(oldpar2)
    dev.off()
  }
}

#########################################################################################
#########################################################################################





## plot viruses

q<-get.centrifuge2(all.list, filter.taxRank=cur.rank, filter.taxon=c("Homo sapiens", "synthetic construct", " phiX174"), 
                   invert.taxon=T, filter.taxonomy=c("Viruses"), absolute=2)
oldpar<-par(mar = par()$mar + c(0,0,0,7))
qq<-barplot(q$m, col=palette()[-1], axisnames=F, main=taxRanks[taxRanks[,1]==cur.rank,2], ylab="proportion of reads", ylim=0:1)
bin<-(qq[2]-qq[1])
q$groups$Ncum<-cumsum(q$groups$N)
f1<-qq[1]-bin/2
bins<-c(f1, q$groups$Ncum*bin+f1)
axis(1, at=bins, labels=rep("",length(bins)))
axis(1, at=qq[q$groups$idx], labels=q$groups$group, las=3, tick=F)
abline(v=bins, col="black", lwd=2)
lines(rep(0,length(bins)))
legend("topright", q$my.legend, cex=1, bty="n", fill=palette()[-1], inset=c(-0.25, -0.05), xpd=T)
par(oldpar)


q<-get.centrifuge(all.list, filter.taxRank="S", filter.taxon=c("Homo sapiens", "synthetic construct", " phiX174"), 
                   invert.taxon=T, filter.taxonomy=c("Viruses"), absolute=2)
q$m$species2<-unlist(lapply(q$my.leg, FUN))
f<-q$m[, lapply(.SD, sum, na.rm=TRUE), by=q$m$species2 ]
q$my.leg<-f$q
q$m<-f[,-1, with=F]
cum.sum<-apply(q$m,1, sum)
my.sort<-order(cum.sum, decreasing=T)
my.legend<-paste0(q$my.leg, " ", format(100*cum.sum/sum(cum.sum, na.rm=T), digits=0, scientific=F),"%")[my.sort]
m<-prop.table(as.matrix(q$m[my.sort,]),2)
oldpar<-par(mar = par()$mar + c(0,0,0,10))
qq<-barplot(m[1:24,], col=palette()[-1], axisnames=F, main=taxRanks[9,2], ylab="proportion of reads", ylim=0:1)
bin<-(qq[2]-qq[1])
q$groups$Ncum<-cumsum(q$groups$N)
f1<-qq[1]-bin/2
bins<-c(f1, q$groups$Ncum*bin+f1)
axis(1, at=bins, labels=rep("",length(bins)))
axis(1, at=qq[q$groups$idx], labels=q$groups$group, las=3, tick=F)
abline(v=bins, col="black", lwd=2)
lines(rep(0,length(bins)))
legend("topright", my.legend[1:24], cex=1, bty="n", fill=palette()[-1], inset=c(-0.51, -0.3), xpd=T)
par(oldpar)

q<-get.centrifuge2(all.list, filter.taxRank="G", filter.taxon=c("Homo sapiens", "synthetic construct", " phiX174"), 
                   invert.taxon=T, filter.taxonomy=c("Viruses"), absolute=2)

x<-all.list[[1]]
a<-x[grep("S_Streptococcus|D_Bacteria", x$taxonomy) ,]$taxonomy
aa<-unlist(strsplit(a, "; "))


for(i in 1:nrow(taxRanks)){
  print(table(aa[grep(paste0(taxRanks[i,1],"_"), aa)]))
}

############################################################################
############################################################################
## check with Anna's stats
ref<-read.table("nb_reads_anna.txt", header=T)
ref<-ref[mixedorder(as.character(ref$Sample)),]
ref$Total<-ref$Total/2

## library sizes
libs<-read.table("library_size.txt", header=F)
libs<-libs[mixedorder(as.character(libs$V1)),]
libs$V2<-libs$V2/4

pdf("Australia_number_of_reads_change.bk.pdf", paper="a4", width=0, height=0)
par(mfrow = c(2,1)) #add room for the rotated labels

plot(ref$Total,libs$V2, xlab="Copenhagen bam files", ylab="Bern Fastq files", main="Copenhagen => Bern")
abline(0,1)

barplot(ref$Total -libs$V2, col=data$group, main="lost reads from the conversion from the bam files to the fastq files")
axis(1, at=m[groups$idx], labels=groups$group, las=3)

plot(libs$V2, data[,1]/2, xlab="Bern Fastq files", ylab="Bern bam files", main="Bern: fastq => bam")
abline(0,1)

barplot(libs$V2 - data[,1]/2, col=data$group, main="lost reads during the mapping")
axis(1, at=m[groups$idx], labels=groups$group, las=3)

dev.off()

plot(ref$Total,data[,1])
abline(0,1)
rownames(data)[abs(ref$Total - data[,1]) > 1e7]
barplot(ref$Total - data[,1])





library_size
