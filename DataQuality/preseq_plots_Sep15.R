setwd("~/Projects/Botocudos/Scripts/")
source("~/Projects/Botocudos/Scripts/get_constructors.R")
source("~/Projects/Botocudos/Scripts/translate_ids.R")
library(parallel)
library(preseqR)
library(ggplot2)
library(scales)
get_dataframes <- function(f, total, Target, extrap=3e8){
  
  # Size of experiments up to extrap sequenced reads
  size_f <- seq(0, round(extrap/total), length.out = 10000)
  
  # Total reads aligning to target regions
  f_total <- sum(f[,1]*f[,2])
  
  # Slopes, equal to "endogenous content", dup included
  slope <- f_total/total

  results <- get_constructor(f, r = 1)
  # Data frames
  
  df <- data.frame(Total_reads = total*size_f, # converts to total reads
                     Expected = tryCatch({sapply(size_f, 
                                                 function(i) 
                                                   results$FUN.bootstrap(i))}, 
                                         error = function(err){
                                           # this = total reads times "endogenous"
                                           return(slope*f_total*size_f/slope)}),
                     Var = sapply(size_f, 
                                  function(i) tryCatch({results$var(i)},
                                                       error = function(err){
                                                         return(0)})))
  reads <- df

  reads$Target <- Target
  
  newList <- list("Reads" = reads, 
                  "results" = results)
  
  return(newList)
}
boto <- read.table("~/Projects/Botocudos/Files/Summaries/Botocudos_summary_2017_09_18.table", header = T)
# Change to expected coverage
###############################################################################

###
# MN00010
f <- read.table("../Files/Histograms/MN00010.txt")
# Total retained reads
total <- 18335055
mn00010 <- list("f"=f, "total"=total, Target = "MN00010")

###
# MN00013
f <- read.table("../Files/Histograms/MN00013.txt")
# Total retained reads
total <- 16620709
mn00013 <- list("f"=f, "total"=total, Target = "MN00013")

###
# MN00016
f <- read.table("../Files/Histograms/MN00016.txt")
# Total retained reads
total <- 10947557
mn00016 <- list("f"=f, "total"=total, Target = "MN00016")

###
# MN00021
# Read files
f <- read.table("../Files/Histograms/MN00021.txt")
# Total retained reads
total <- 9551854
mn00021 <- list("f"=f, "total"=total, Target = "MN00021")

###
# MN00022
# Read files
f <- read.table("../Files/Histograms/MN00022.txt")
# Total retained reads
total <- 4324624
mn00022 <- list("f"=f, "total"=total, Target = "MN00022")

###
# MN00023
# Read files
f <- read.table("../Files/Histograms/MN00023.txt")
# Total retained reads
total <- 2631186
mn00023 <- list("f"=f, "total"=total, Target = "MN00023")

###
# MN0003
# Read files
f <- read.table("../Files/Histograms/MN0003.txt")
# Total retained reads
total <- 7652338
mn0003 <- list("f"=f, "total"=total, Target = "MN0003")

###
# MN00039
f <- read.table("../Files/Histograms/MN00039.txt")
# Total retained reads
total <- 27452692
mn00039 <- list("f"=f, "total"=total, Target = "MN00039")

###
# MN00045
f <- read.table("../Files/Histograms/MN00045.txt")
# Total retained reads
total <- 9832730
mn00045 <- list("f"=f, "total"=total, Target = "MN00045")

###
# MN00019
f <- read.table("../Files/Histograms/MN00019.txt")
# Total retained reads
# lane1.1:2613840
# lane1.2:2618040
# lane1.3:2622055
# lane1.4:2612130
# lane1.5:817122
# Nuclear: 8416
# 

total <- 2613840 + 2618040 + 2622055 + 2612130 + 817122
mn00019 <- list("f"=f, "total"=total, Target = "MN00019")

###
# MN00056
# Read files
f <- read.table("../Files/Histograms/MN00056.txt")
# Total retained reads
total <- 2341638
mn00056 <- list("f"=f, "total"=total, Target = "MN00056")

###
# MN00064
f <- read.table("../Files/Histograms/MN00064.txt")
# Total retained reads
total <- 14325117
mn00064 <- list("f"=f, "total"=total, Target = "MN00064")

###
# MN00066
f <- read.table("../Files/Histograms/MN00066.txt")
# Total retained reads
total <- 13287914
mn00066 <- list("f"=f, "total"=total, Target = "MN00066")

###
# MN00067
# Read files
f <- read.table("../Files/Histograms/MN00067.txt")
# Total retained reads
total <- 5993189
mn00067 <- list("f"=f, "total"=total, Target = "MN00067")

### 
# MN00068
f <- read.table("../Files/Histograms/MN00068.txt")
# Total retained reads
total <- 9680238
mn00068 <- list("f"=f, "total"=total, Target = "MN00068")

###
# MN00069
f <- read.table("../Files/Histograms/MN00069.txt")
# Total retained reads
total <- 7757426
mn00069 <- list("f"=f, "total"=total, Target = "MN00069")

###
# MN0009
f <- read.table("../Files/Histograms/MN0009.txt")
# Total retained reads
total <- 7064993
mn0009 <- list("f"=f, "total"=total, Target = "MN0009")

###
# MN0118
f <- read.table("../Files/Histograms/MN00118.txt")
# Total retained reads
total <- 11622382
mn00118 <- list("f"=f, "total"=total, Target = "MN00118")

###
# MN0119
f <- read.table("../Files/Histograms/MN00119.txt")
# Total retained reads
total <- 7523029
mn00119 <- list("f"=f, "total"=total, Target = "MN00119")


###
# MN00316
f <- read.table("../Files/Histograms/MN00316.txt")
# Total retained reads
total <- 10732101
mn00316 <- list("f"=f, "total"=total, Target = "MN00316")

###
# MN00346
f <- read.table("../Files/Histograms/MN00346.txt")
# Total retained reads
total <- 13755034
mn00346 <- list("f"=f, "total"=total, Target = "MN00346")

###
# MN1701
f <- read.table("../Files/Histograms/MN01701.txt")
# Total retained reads
total <- 13349357
mn01701 <- list("f"=f, "total"=total, Target = "MN01701")

###
# MN1943
# Read files
f <- read.table("../Files/Histograms/MN1943.txt")
# Total retained reads
total <- 14291106
mn1943 <- list("f"=f, "total"=total, Target = "MN1943")


#########
#.........................................................................
#

my_hist <- list(MN00010 = mn00010, MN00013 = mn00013,
                MN00016 = mn00016, MN00021 = mn00021,
                MN00022 = mn00022, MN00023 = mn00023,
                MN0003 = mn0003, MN00039 = mn00039,
                MN00045 = mn00045, MN00019 = mn00019,
                MN00056 = mn00056, MN00064 = mn00064,
                MN00066 = mn00066, MN00067 = mn00067,
                MN00068 = mn00068, MN00069 = mn00069,
                MN0009 = mn0009, MN00118 = mn00118,
                MN00119 = mn00119, MN00316 = mn00316,
                MN00346 = mn00346, MN017101 = mn01701,
                MN1943 = mn1943)


# Calculate the number of cores
no_cores <- detectCores() - 1

# Initiate cluster
cl <- makeCluster(no_cores)

clusterExport(cl, list("get_dataframes", "my_hist", "get_constructor"))
clusterEvalQ(cl, library(preseqR))

semi_final <- parLapply(cl, my_hist, function(target){
  get_dataframes(f = target$f,
                 total = target$total, Target = target$Target)
})
stopCluster(cl)
save(semi_final, file = "~/Projects/Botocudos/Files/PreSeq/output_preseq_Sep16_all.Rda")
system("say -v Thomas Il est finit")

#------------------------------------------------------------------------------
# Plots

# Colors
colors <- colorMe(final$STM2$Reads$Condition)
Cond_labels <- c("Pre-capture", "Y", "WGC", "WGC+Y")
##### Complexity curves
plots <- list()
#######
## Save to pdf
plots <- list()
for(n in names(semi_final)){
  p <- ggplot(semi_final[[n]]$Reads, 
         aes(x = Total_reads, 
             y = Expected,
             ymin = Expected - sqrt(Var),
             ymax = Expected + sqrt(Var))) +
    geom_line() +
    geom_ribbon(alpha = 0.6, fill = "chartreuse3")  +
    geom_vline(xintercept = my_hist[[n]]$total, linetype = "dashed", color = "gray") +
    labs(title = paste(n, " (", 
                       percent(signif(boto$hits_unique_frac.endogenous.[
                         boto$Library == n], 2)), ")", sep = ""), 
         x = "Total sequenced reads", y = "Expected unique reads") +
    theme(legend.position = "none")
  plots[[length(plots)+1]] <- p
}

pdf("~/Projects/Botocudos/Plots/Preseq_expected_reads.pdf", width = 15, height = 15)
png("~/Projects/Botocudos/Plots/Preseq_expected_reads.png", 
    width = 12, height = 30, units = "in", res = 250)
plot_grid(plotlist = plots, ncol = 3)
dev.off()

############
# Genome coverage
plots <- list()
for(n in names(semi_final)){
  if(n == "MN017101"){next}
  m <- mn2ma(n)
  f <- read.table(paste("~/Projects/Botocudos/Files/PreSeq/GC/", 
                       m,".hg19_nuc.realigned.bam.mr_gc_extrap.txt", sep = ""),
                  header = T)
  endo <- boto$hits_unique_frac.endogenous.[boto$Library == n]
  l <- boto$hits_length.nuclear.[boto$Library == n]
  p <- ggplot(f,
              aes(x = TOTAL_BASES/l*endo, 
                  y = (EXPECTED_COVERED_BASES)/(1e+9),
                  ymin = (EXPECTED_COVERED_BASES - LOWER_95.CI)/(1e+9),
                  ymax = (EXPECTED_COVERED_BASES + LOWER_95.CI)/(1e+9))) +
    geom_line() +
    geom_ribbon(alpha = 0.6, fill = "salmon")  +
    coord_cartesian(xlim = c(0, 1e+7)) +
    #geom_vline(xintercept = my_hist[[n]]$total, linetype = "dashed", color = "gray") +
    labs(title = paste(n, " (", 
                       percent(signif(boto$hits_unique_frac.endogenous.[
                         boto$Library == n], 2)), ")", sep = ""), x = "Total sequenced reads", y = "Fraction of the genome") +
    theme(legend.position = "none")
  plots[[length(plots)+1]] <- p
  
}

pdf("~/Projects/Botocudos/Plots/Preseq_gc.pdf", width = 10, height = 20)
png("~/Projects/Botocudos/Plots/Preseq_gc.png", 
    width = 12, height = 30, units = "in", res = 250)

plot_grid(plotlist = plots, ncol = 3)
dev.off()

###########################################################
#### January 31st
boto2 <- read.table("~/Projects/Botocudos/Files/Summaries/Botocudos_summary_2017_12_06.table",
                    header = T)
my_hist <- list()

individuals <- unique(boto$Target)
i = 1
for(ind in individuals){
  f <- read.table(paste("~/Projects/Botocudos/Files/Histograms/2017_12_06/Merged/",
                        ind, ".bam_duphist.txt", sep = ""))
  total <- unique(boto$seq_reads_se[boto$Target == ind] + boto2$seq_reads_se[boto2$Target == ind])
  x <- list("f"=f, "total"=total, Target = ma2mn(ind))
  my_hist[[ind]] <- x
  i <- i + 1
}

# Calculate the number of cores
no_cores <- detectCores() - 1

# Initiate cluster
cl <- makeCluster(no_cores)

clusterExport(cl, list("get_dataframes", "my_hist", "get_constructor"))
clusterEvalQ(cl, library(preseqR))

semi_final <- parLapply(cl, my_hist, function(target){
  get_dataframes(f = target$f,
                 total = target$total, Target = target$Target, extrap = 10e8)
})
stopCluster(cl)
save(semi_final, file = "~/Projects/Botocudos/Files/PreSeq/output_preseq_Jan31_all.Rda")
system("say -v Thomas Il est finit")
