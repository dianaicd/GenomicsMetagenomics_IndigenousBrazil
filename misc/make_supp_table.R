# Script to build a supplementary table with statistics related to mapping
#------------------------------------------------------------------------------#
# Load libraries and sources
#------------------------------------------------------------------------------#
library(plyr)
library(ggplot2)
library(scales)
library(cowplot)
library(ggridges)
source("~/Projects/Botocudos/Scripts/Plotting/length_distribution_plot.R")
source("~/Projects/Botocudos/Scripts/Plotting/mapDamage_plot.R")
source("~/Projects/Botocudos/Scripts/Complexity/preseq_functions.R")

#------------------------------------------------------------------------------#
# Parse arguments
args <- commandArgs(T)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1

set_def_arg <- function(arg, def, given){
  
}

## Arg1 default
if(is.null(args$arg1)) {
  ## do something
}

## Arg2 default
if(is.null(args$arg2)) {
  ## do something
}

#------------------------------------------------------------------------------#
# Order in which DNA libraries will be displayed
lib_order <- c("All", "mtCapture", "L1", "L2", "L3U")
colors_per_library <- c("gray40", "#c4bf33", "#b08699","#e0b81f","#0d75ff")
names(colors_per_library) <- lib_order
paste_path <- function(end, path = "~/Projects/Botocudos/Plots/"){
  return(paste(path, end, sep = ""))
}

#------------------------------------------------------------------------------#
# Read in summary files
myDir <- "~/Projects/Botocudos/Files/Summaries/2019_07_26/"

files <- list.files(path = myDir, pattern = ".new.summary$")

boto <- data.frame()

for(file in files){
  tmp_summary <- read.table(paste(myDir, file, sep =""), header = T, sep = "\t")
  boto <- rbind(boto, tmp_summary)
}

boto$library <- sub("\\*", "All", boto$library)
boto$library <- sub("S1", "L1", boto$library)

boto$library <- factor(boto$library, levels = lib_order, ordered = T)
#------------------------------------------------------------------------------#
# Merge isotopic data if available
backbone <- read.csv("~/Projects/Botocudos/Files/Tables/Backbone_Table1.csv")

colnames(backbone) <- c("sample", "CalibratedDate", "Area", "State") 

boto <- join(boto, backbone, by = "sample")

boto <- boto[order(boto$hits_unique_frac_endogenous, boto$library),]
boto$sample <- factor(boto$sample, levels = unique(boto$sample), ordered = T)
#------------------------------------------------------------------------------#
# Merge contamination data if available
#------------------------------------------------------------------------------#
## MT
prefix <- "~/Projects/Botocudos/Files/Contamination/2019_08_06/"
samples <-unique(boto$sample)
rmTrans <- c("all", "rmTrans")


rmBurnIn_contammmix <- function(path, sample, rmTrans, size = 100, lib = "All",
                                 burnIn = 0.1){
  if(lib == "All"){
      myRda <- paste(path, sample, "/",sample, "_", rmTrans, ".Rdata", sep = "")
  }else{
    myRda <- paste(path, sample, "/", lib, "_", rmTrans, ".Rdata", sep = "")
  }
  print(myRda)
  load(myRda)
  notBurnIn <- as.integer(burnIn*size)
  tiny <-
    #remove burnin
    c(res$chains[[1]][size - ((size - notBurnIn):0),2])
                 
  result <- data.frame(estimate = tiny, damage = rmTrans, lib = lib,
                       sample = sample)
  return(result)
}

if(load_contammix){
  load("~/Projects/Botocudos/Files/Contamination/2019_08_06/contamination.Rdata")
}else{
  contammix <- data.frame()
  
  for(s in samples){
    libs <- unique(boto$library[boto$sample == s])
    for(l in libs){
      for(r in rmTrans){
        tmp_contam <- rmBurnIn_contammmix(prefix, s, r, 100000, l)
        contammix <- rbind(contammix, tmp_contam)
      }
    }
  }
}

# Extract mitochondrial coverage and number of reads used in estimation
path <- "~/Projects/Botocudos/Files/Contamination/2019_08_06/"
samples <- unique(boto$sample)
rmTrans <- c("all", "rmTrans")
annot_mito <- data.frame()

for(s in samples){
  libs <- unique(boto$library[boto$sample == s])
  for(lib in libs){
    for(r in rmTrans){
      if(lib == "All"){
        myRda <- paste(path, s, "/",s, "_", r, ".Rdata", sep = "")
      }else{
        myRda <- paste(path, s, "/", lib, "_", r, ".Rdata", sep = "")
      }
      #print(myRda)
      load(myRda)
      nReads <- dim(res$mnMatrix)[1]
      coverage <- boto$hits_coverage_mitochondrial[boto$sample == s & boto$library == lib]
      annot_tmp <- data.frame(nReads = nReads, coverage = coverage, damage = r, lib = lib,
                              sample = s)
      
      annot_mito <- rbind(annot_mito, annot_tmp)
    }
  }
}

MAP <- read.csv("~/Projects/Botocudos/Files/Contamination/2019_08_06/estimates.txt", header = F)

colnames(MAP) <- c("sample", "lib", "damage","map", "low", "high")
MAP$lib <- sub("S1", "L1", MAP$lib)
annot_mito <- join(annot_mito, MAP, by = c("sample", "lib", "damage"))
#------------------------------------------------------------------------------#
## X
minDepth <- seq(1,3)
path <- "~/Projects/Botocudos/Files/Contamination/2019_08_06/X/"
mySamples <- unique(boto$sample)
libs <- unique(boto$library)

s <- mySamples[1]
l <- libs[2]
md <- minDepth[3]

empty_df <- function(sample, lib){
  myDF <- data.frame(method = c("One-cns", "Two-cns"),
                     estimate = c(NA, NA), lb = c(NA, NA), ub = c(NA, NA),
                     err = c(NA, NA), nSites = c(NA, NA))
  return(myDF)
}

x_cont <- data.frame()
for(s in mySamples){
  for(l in libs){
    for(md in minDepth){
      if(l == "All"){
        tmp_res <- tryCatch(tmp_res <- read.table(paste(path, s, "/", s, 
                                                        "_md", md, ".result", sep = "")),
                            error = function(e) empty_df(s, l),
                            finally = function() return(tmp_res)
        )
      }else{
        tmp_res <- tryCatch(tmp_res <- read.table(paste(path, s, "/", l, 
                                                        "/", l, "_md", md, ".result", sep = "")),
                            error =  function(e) empty_df(s, l),
                            finally = function() return(tmp_res))
      }
      colnames(tmp_res) <- c("method", "estimate", "lb", "ub", "err", "nSites")
      tmp_res$sample <- s
      tmp_res$library <- l
      tmp_res$depth <- md
      x_cont <- rbind(x_cont, tmp_res)
    }
  }
}

x_cont <- join(x_cont, boto[, c("sample", "library", "sex")], by = c("sample", "library"))
x_cont <- x_cont[!is.na(x_cont$sex),]

males <- x_cont[x_cont$sex %in% c("XY", "consistent with XY but not XX"),]

males <- males[males$method == "Two-cns",]
males <- males[order(males$sample, males$library),]

# write.table(males, "~/Projects/Botocudos/Files/Contamination/2019_08_06/X/estimates_X_males.txt", 
# sep = "\t", quote = F, col.names = T, row.names = F)

males_fmt <- males

males_fmt$estimate <- percent(males_fmt$estimate)
males_fmt$lb <- percent(males_fmt$lb)
males_fmt$ub <- percent(males_fmt$ub)

#males_fmt
males_dp3 <- males[males$depth == 3,]
colnames(annot_mito) <- c("mito_num_Reads", "mito_coverage", "mito_polymorphic_sites",
                          "library", "sample", "mito_estimate", "mito_low_est",
                          "mito_upper_est")
colnames(males_dp3) <- c("X_method", "X_estimate", "X_low_est", "X_upper_est",
                         "X_error", "X_num_sites", "sample", "library", "X_min_depth",
                         "sex")
estimates <- join(males_dp3, annot_mito[annot_mito$mito_polymorphic_sites == "all",],
                  by = c("sample", "library"))

estimates <- estimates[,c("sample", "library",
                          "mito_estimate", "mito_low_est", "mito_upper_est",
                          "X_estimate", "X_low_est", "X_upper_est",
                          "mito_num_Reads", "mito_coverage",
                          "X_num_sites", "X_min_depth",
                          "sex")]

estimates$mito_estimate <- 1 - estimates$mito_estimate
estimates$mito_low_est <- 1 - estimates$mito_low_est
estimates$mito_upper_est <- 1 - estimates$mito_upper_est
tmp_cont <- estimates$mito_low_est
estimates$mito_low_est <- estimates$mito_upper_est
estimates$mito_upper_est <- tmp_cont

for(myColumn in c("mito_estimate", "mito_low_est", "mito_upper_est",
                  "X_estimate", "X_low_est", "X_upper_est",
                  "mito_coverage") ){
  estimates[, myColumn] <- signif(estimates[,myColumn], 3)
}

for(myColumn in c("mito_estimate", "mito_low_est", "mito_upper_est",
                  "X_estimate", "X_low_est", "X_upper_est"
) ){
  estimates[, myColumn] <- percent(estimates[,myColumn], 2)
}

#------------------------------------------------------------------------------#
# Do the excel shit
mito <- MAP
colnames(mito) <- c("sample", "library","damage", "cont_MT", "cont_MT_upper", "cont_MT_lower")
mito[,c("cont_MT", "cont_MT_upper", "cont_MT_lower")] <- 1 - mito[,c("cont_MT", "cont_MT_upper", "cont_MT_lower")]
mito <- mito[mito$damage == "all",]
mito$damage <- NULL

X <- males_dp3
colnames(X) <- c("X_method", "cont_X","cont_X_lower", "cont_X_upper", "X_err","X_sites_dp3", "sample", "library", "depth", "sex")
X <- X[,c("sample", "library","cont_X", "cont_X_lower", "cont_X_upper")]


mito[,c("cont_MT", "cont_MT_upper", "cont_MT_lower")] <- round( mito[,c("cont_MT", "cont_MT_upper", "cont_MT_lower")], 2)
X[, c("cont_X", "cont_X_lower", "cont_X_upper")] <- round(X[, c("cont_X", "cont_X_lower", "cont_X_upper")], 2)

boto_complete <- join(boto, mito, by = c("sample", "library"))
boto_complete <- join(boto_complete, X, by = c("sample", "library"))



boto_complete$cont_MT <- paste(percent(boto_complete$cont_MT), 
                               " (", percent(boto_complete$cont_MT_lower), 
                               ", ", percent(boto_complete$cont_MT_upper), ")",
                               sep = "")

boto_complete$cont_X <- paste(percent(boto_complete$cont_X), 
                              " (", percent(boto_complete$cont_X_lower), 
                              ", ", percent(boto_complete$cont_X_upper), ")",
                              sep = "")
boto_complete$seq_retained_frac <- 1 - boto_complete$seq_trash_se_frac

to_round <- as.data.frame(matrix(data=c("Ry", 3,
                                        "Ry_confint", 3,
                                        "seq_trash_se_frac", 2,
                                        "seq_retained_length", 1,
                                        "hits_raw_frac_endogenous", 4,
                                        "hits_unique_frac_endogenous", 4, 
                                        "hits_coverage_endogenous", 3,
                                        "hits_length_endogenous", 1,
                                        "hits_unique_frac_mitochondrial", 6,
                                        "hits_coverage_mitochondrial",1,
                                        "hits_length_mitochondrial", 1,
                                        "hits_length_nuclear",1,
                                        "hits_coverage_nuclear", 3),
                                 ncol = 2, byrow = T), stringsAsFactors = F)
to_round$V2 <- as.integer(to_round$V2)

#to_round_digits <- c(3,3,2,1,4,4,3,1,6,1,1,3)
for(i in 1:nrow(to_round)){
  boto_complete[,to_round[i,1]] <- round(boto_complete[, to_round[i,1]], 
                                         to_round[i,2])
  
}
#------------------------------------------------------------------------------#
# Compute other numbers to report in the table
boto_complete$hits_clonality_endogenous_frac <- percent(boto_complete$hits_clonality_endogenous/boto_complete$hits_raw_endogenous)

boto_complete$hits_raw_frac_endogenous <- percent(boto_complete$hits_raw_frac_endogenous, accuracy = 0.01)
boto_complete$hits_unique_frac_endogenous <- percent(boto_complete$hits_unique_frac_endogenous, accuracy = 0.01)

#boto_complete$hits_raw_frac_nuclear <- percent(boto_complete$hits_raw_frac_nuclear, accuracy = 0.01)
boto_complete$hits_unique_frac_nuclear <- percent(boto_complete$hits_unique_frac_nuclear, accuracy = 0.01)

boto_complete$hits_unique_frac_mitochondrial <- percent(boto_complete$hits_unique_frac_mitochondrial, accuracy = 0.0001)
boto_complete$seq_trash_se_frac <- percent(boto_complete$seq_trash_se_frac)
boto_complete$seq_retained_frac <- percent(boto_complete$seq_retained_frac)

boto_complete$Ry <- paste(boto_complete$Ry, " (", boto_complete$Ry-boto_complete$Ry_confint, ", ",
                          boto_complete$Ry+boto_complete$Ry_confint, ")", sep = "")
#------------------------------------------------------------------------------#
# Add thousands separator
need_comma <- which(unlist(sapply(1:ncol(boto_complete), 
                                  function(x) class(boto_complete[,x])) == "integer"))

boto_complete[,need_comma] <- sapply(need_comma, 
                                     function(x) prettyNum(boto_complete[,x], big.mark = ","))

myColNames <- c("sample",                         "Sample", "National Museum of Brazil (Museu Nacional, MN) identifier",
                # "target",                         "Sample2",
                "Area",                           "Region", "Sample origin (region) according to MN's archives",
                "State",                          "State", "Sample origin (state) according to MN's archives", 
                "CalibratedDate",                 "Calibrated date (AD)", "OxCal estimated mean and standard deviation for calibrated dates",
                
                "library",                        "Library", "DNA library identifier",
                "lib_type",                       "library_type", "Library type",
                
                "seq_reads_se",                   "sequenced_reads", "Number of sequenced reads",
                #"seq_trash_se",                   "trashed_reads_trimming", 
                #"seq_trash_se_frac",              "fraction_trashed_reads_trimming",
                "seq_retained_reads",             "retained_reads_trimming", "Number of reads retained after trimming",
                "seq_retained_frac",              "percent_retained_reads_trimming", "Percentage of sequenced reads retained after trimming",
                #"seq_retained_nts",               "retained_nucleotides_trimming",
                "seq_retained_length",            "length_retained_reads", "Average read length for retained reads after trimming",
                
                "hits_clonality_endogenous",      "duplicated_reads", "Number of reads flagged as PCR duplicates by picardtools MarkDuplicates",
                "hits_clonality_endogenous_frac", "percent_duplicated_reads", "Percentage of reads flagged as PCR duplicates by picardtools MarkDuplicates",
                
                "hits_unique_endogenous",         "unique_reads_endogenous", "Number of unique reads mapped to the human genome (build 37.1)",
                "hits_unique_frac_endogenous",    "percent_unique_reads_endogenous", "Percentage of retained reads that were uniquely mapped to the human genome",
                "hits_length_endogenous",         "length_reads_endogenous", "Average length of uniquely mapped reads",
                "hits_coverage_endogenous",       "DoC_endogenous", "Depth of coverage on human genome (number of unique bases divided by the length of the human genome reference build 37.1)",
                
                #"hits_unique_nuclear",            "unique_reads_nuclear",
                #"hits_coverage_nuclear",          "DoC_nuclear",
                "hits_unique_mitochondrial",      "unique_reads_MT", "Number of unique reads mapped to the mitochondrial genome",
                "hits_coverage_mitochondrial",    "DoC_MT", "Depth of coverage on mitochondrial genome (number of unique bases divided by the length of the mitochondrial genome, build 37.1)",
                #"hits_raw_endogenous",            "reads_endogenous_raw",
                #"hits_raw_frac_endogenous",       "fraction_reads_raw_endogenous",
                
                
                #"hits_unique_frac_nuclear",       "percent_unique_reads_nuclear",
                #"hits_length_nuclear",            "length_reads_nuclear",
                
                #"hits_unique_frac_mitochondrial", "fraction_unique_reads_mitochondrial",
                "hits_length_mitochondrial",      "length_reads_mitochondrial", "Average length of uniquely mapped reads on mitochondrial genome",
                
                "cont_MT",                        "contamination_MT", "Maximum a posteriori estimate and 95% credible interval for contamination based on the mitochondrial genome (contamMix)",
                #"cont_MT_upper",                  "contamination_mitochondrial_upper",
                #"cont_MT_lower",                  "contamination_mitochondrial_lower",
                "cont_X",                         "contamination_X", "Maximum likelihood estimate and 95% confidence interval for contamination based on X-chromosome data (contaminationX)",
                #"cont_X_lower",                   "contamination_X_lower",
                #"cont_X_upper",                   "contamination_X_upper",
                
                "Ry",                             "Ry", "Ratio of the reads mapping to the Y-chromosome and the reads mapping to sexual chromosomes (X and Y) and 95% confidence interval",
                #"Ry_confint",                     "Ry_confint",
                "sex",                            "Sex", "Molecular sex determined by the Ry ratio"
)

change_name <- function(name){
  new_name <- myColNames[which(myColNames == name)[1] + 1]
  ifelse(is.na(new_name), return(name), return(new_name))
}


colnames(boto_complete) <- sapply(colnames(boto_complete), function(name) change_name(name))

boto_complete <- boto_complete[,myColNames[seq(2, length(myColNames), 3)]]

boto_complete$Library <- factor(boto_complete$Library, levels = c("mtCapture",
                                                                  "L1",
                                                                  "L2",
                                                                  "L3U",
                                                                  "All"),
                                ordered = T)

boto_final <- subset(boto_complete, select = myColNames[seq(2, length(myColNames), 3)])
#------------------------------------------------------------------------------#

# boto_complete <- boto_complete[order( boto_complete$Sample, boto_complete$Library,boto_complete$percent_unique_reads_endogenous ),]
boto_final <- boto_final[order(boto_final$Sample, boto_final$Library, boto_final$percent_unique_reads_endogenous),]
write.table(boto_final, "~/Projects/Botocudos/Files/Summaries/2019_07_26/Sup_Table_1.txt",
            col.names = T, row.names = F, sep = "\t", quote = F)
#------------------------------------------------------------------------------#
require(xlsx)
wb_name <- "tabla1.xlsx"
sh_name <- "SI_Table1"
# define workbook and sheet
wb <- createWorkbook(type="xlsx")
sheet <- createSheet(wb, sheetName = sh_name)

# Define row style to be used
mySamples <- data.frame(sample = levels(boto_final$Sample), 
                        color = rep(c("blue", "white"), 
                                    round(length(unique(boto_final$Sample))/2)))
white_row <- CellStyle(wb) 
blue_row <- CellStyle(wb) + Fill(foregroundColor = "lightblue")

# Function to add a row with specified cell style
# xlsx.addRow <- function(sheet, rowIndex, data, rowStyle, sample){
#   rows <- createRow(sheet, rowIndex = rowIndex)
#   sheetRow <- createCell(rows, colIndex=1)
#   sapply(data, function(x) setCellValue(sheetRow[[1,1]], x))
#   cells <- getCells(getRows(sheet, rowIndex = rowIndex))
#   
#   if(mySamples$color[mySamples$sample == sample] == "blue"){
#      lapply(names(cells), function(ii) setCellStyle(sheetRow[[ii]], blue_row))
#   }else{
#     #setCellStyle(sheetRow[[1,1]], white_row)
#   }
# }

# for(i in 1:nrow(boto_final)){
#   xlsx.addRow(sheet = sheet, rowIndex = i, data = boto_final[i,],
#               sample = as.character(boto_final$Sample[i]))
# }

# Add table + column name explanation
addDataFrame(boto_final, sheet, startRow = 1, startColumn = 1, row.names = F)

# Add color
rows <- getRows(sheet, rowIndex = 2:(nrow(boto_final)+1))
cells <- getCells(rows, colIndex = 1:ncol(boto_final))
values <- lapply(cells, getCellValue)

highlight_blue <- c()
for(i in names(values)){
  if(values[i] %in% mySamples$sample){
    currentSample <- values[i]
  }
  if(mySamples$color[mySamples$sample == currentSample] == "blue"){
    highlight_blue <- c(highlight_blue, i)
  }
}
lapply(names(cells[highlight_blue]), function(ii) setCellStyle(cells[[ii]], blue_row))



myLegend <- as.data.frame(matrix(myColNames, ncol = 3, byrow = T) )[,2:3]; #print(myLegend)
addDataFrame(myLegend, sheet, startRow = nrow(boto_final) + 3,
             startColumn = 1, row.names = F, col.names = F)

autoSizeColumn(sheet, colIndex = 1:ncol(boto_final))
saveWorkbook(wb, wb_name)
