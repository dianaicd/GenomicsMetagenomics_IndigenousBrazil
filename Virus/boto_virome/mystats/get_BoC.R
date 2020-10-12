#!/usr/bin/env Rscript

# August 14 2020.
# Script to calculate the breadth of coverage given a bedtools coverage table.

# August 26 2020.
# Get depth of coverage, and standard deviations only considering sites with at least one read 
# (SD BoC) & one considering all sites (SD DoC) 

# October 9 2020
# Get the effective DoC (DoC of sites covered at least once).
# What we called SD_BoC it's actually the SD of the effective DoC.

# Actions:
# Read coverage table.
# Calculate breadth of coverage, depth of coverage, effective depth of coverage & standard deviations of DoC & eff. DoC
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
SD_DoC <- bedtools_coverage[, sd(V3)]
eff_DoC <- bedtools_coverage[, sum(V3)]/bedtools_coverage[V3>0, .N]
# When no read maps
eff_DoC <- ifelse(is.na(eff_DoC), 0, eff_DoC)
SD_eff_DoC <- bedtools_coverage[, sd(V3>0)]

final_table <- data.table(BoC = BoC, 
                          DoC = DoC,
                          SD_DoC = SD_DoC,
                          eff_DoC = eff_DoC,
                          SD_eff_DoC = SD_eff_DoC
                          )

fwrite(final_table, file = output_file, sep = "\t", col.names = F)