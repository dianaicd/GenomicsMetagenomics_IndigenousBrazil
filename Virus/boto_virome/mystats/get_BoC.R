#!/usr/bin/env Rscript

# August 14 2020.
# Script to calculate the breadth of coverage given a bedtools coverage table.

# Actions:
# Read coverage table.
# Calculate breadth of coverage & save it to a file.

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
              help = "Output file")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

cov_table <- opt$cov_table
output_file <- opt$output

bedtools_coverage <- fread(input = cov_table)

BoC <- bedtools_coverage[V3>0, .N]/nrow(bedtools_coverage)

write(BoC, file = output_file)