#!/usr/bin/env Rscript

# August 14 2020.
# Script to calculate the breadth of coverage given a bedtools coverage table.

# August 26 2020.
# Get depth of coverage, and standard deviations only considering sites with at least one read 
# (SD BoC) & one considering all sites (SD DoC) 

# Actions:
# Read coverage table.
# Calculate breadth of coverage, depth of coverage & standard deviations of BoC & DoC
# Save to file.

# Running example:
# Rscript --vanilla get_BoC.R --cov_table {input} --output {output}

#### Libraries ####
library(data.table)
library(optparse)

#### Arguments ####
option_list = list(
  make_option(c("-c", "--cov_table"), type = "character", default = NULL, 
              help = "bedtools coverage table", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = NULL, 
              help = "output file name", metavar = "character")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

cov_table <- opt$cov_table
output_file <- opt$output

bedtools_coverage <- fread(input = cov_table)

BoC <- bedtools_coverage[V3>0, .N]/nrow(bedtools_coverage)
DoC <- bedtools_coverage[, sum(V3)]/nrow(bedtools_coverage)
SD_BoC <- bedtools_coverage[, sd(V3>0)]
SD_DoC <- bedtools_coverage[, sd(V3)]

final_table <- data.table(BoC = BoC, 
                          DoC = DoC, 
                          SD_BoC = SD_BoC, 
                          SD_DoC = SD_DoC)

fwrite(final_table, file = output_file, sep = "\t", col.names = F)