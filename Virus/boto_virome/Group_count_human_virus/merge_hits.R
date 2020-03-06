#!/usr/bin/env Rscript

# February 25 2020
# Script to merge all *hits_by_virus.csv tables from "count_hits_by_virus.R" into a single table.
# Tables to merge must have the suffix "hits_by_virus.csv".
# All *hits_by_virus.csv tables must be in the same folder as the current script.

# February 28 2020: now the suffix of the files must be given to the script

# Actions:
# 1) Read all tables
# 2) Merge the tables

# Running example:
# Rscript --vanilla merge_hits.R -s _hits_by_virus.csv -o all_samples_hits_by_virus.csv

# Print help:
# Rscript --vanilla merge_hits.R -h

#### Libraries: ####
library(optparse)
library(data.table)
library(dplyr)
library(stringr)

#### Arguments list ####
option_list = list(
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Name of the output file", metavar = "character"),
  make_option(c("-s", "--suffix"), type = "character", default = NULL, 
              help = "files' suffix", metavar = "character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

output <- opt$output
suffix <- opt$suffix

#### Open and save the tables ####
hits_by_virus <- list.files(pattern = paste0("*", suffix))
names_tables <- str_remove(hits_by_virus, suffix)
names_tables <- paste0(rep("_", length(names_tables)), names_tables)
hbv_tables <- lapply(hits_by_virus, fread)

#### Merge all tables ####

counter <- 0
all_samples_hits_by_virus <- Reduce(function(...){
  counter <<- counter+1
  merge(..., 
        all = TRUE, 
        by = c("v_names", "tax_id"), 
        suffixes = c(names_tables[counter], names_tables[counter+1]))
  }, 
  hbv_tables)

# Output the table
write.csv(all_samples_hits_by_virus, file = output, quote = F, row.names = F, na = "0")

