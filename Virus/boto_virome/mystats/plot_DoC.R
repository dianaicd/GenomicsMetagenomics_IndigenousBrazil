#!/usr/bin/env Rscript

# September 18 2020.
# Script to create a depth of coverage plot. It recieves as input a coverage table (the one obtained with
# mystats Snakefile)

# Actions:
# Read coverage table
# Plot the depth of coverage & save it to png

# Running example:
# Rscript --vanilla plot_DoC.R --cov_table Bot08.SARS_CoV2.mapq17.coverage.tsv --plot_name
# Bot08.SARS_CoV2.mapq17.coverage.png


#### Libraries ####
library(data.table)
library(optparse)

#### Arguments ####

option_list = list(
  make_option(c("-c", "--cov_table"), type = "character", default = NULL, 
              help = "path to coverage table", metavar = "character"),
  make_option(c("-p", "--plot_name"), type = "character", default = NULL, 
              help = "output plot name", metavar = "character")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

cov_table <- opt$cov_table
plot_name <- opt$plot_name

#### Functions ####
# Plotting function
plot_cov <- function(coverage, plot_name){
  png(plot_name, width = 7, height = 10, units = "in", res = 300)
  
  barplot(height = coverage$V3, 
          main = "Depth of coverage", 
          ylab = "Number of reads",
          col = "#d60036",
          border = NA,
          axes = T,
          space = 0)
  
  xtick <- ceiling(seq(from = 1, to = nrow(coverage), along.with = c(1:5)))
  axis(side = 1, at = xtick, labels = F, tck = -0.005)
  text(x = xtick, par("usr")[3], labels = xtick, pos = 1, xpd = T, cex = 0.45, offset = 0.19, srt = -45)

}

#### Main ####

coverage_data <- fread(cov_table, drop = 1)

plot_cov(coverage = coverage_data, plot_name = plot_name)



