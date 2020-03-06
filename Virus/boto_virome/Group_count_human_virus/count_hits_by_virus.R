#!/usr/bin/env Rscript

# February 20 2020
# Script to get the names of the virus using the taxid & count the hits by virus.

# Actions:
# 1) Read output of python script "select_hvirus_hits.py" (after perl substitution).
# 2) Get virus names with taxize.
# 3) Group & count hits according to virus.

# Running example:
# Rscript --vanilla count_hits_by_virus.R

# Print help:
# Rscript --vanilla count_hits_by_virus.R -h

#### Libraries: ####
library(optparse)
library(data.table)
library(taxonomizr)
library(dplyr)

#library("taxize", lib.loc = "~/R/x86_64-koji-linux-gnu-library/3.5")

# To avoid problems with taxize
# httr::set_config(httr::config(http_version = 0))

#### Arguments list ####
option_list = list(
  make_option(c("-f", "--csv_file"), type = "character", default = NULL, 
              help = "Input csv file name", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Name of the output file", metavar = "character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

csv_file <- opt$csv_file
output <- opt$output

#### Get virus names ####

# Open the output of python script "select_hvirus_hits.py"
hvs_hits <- fread(file = csv_file, sep = ",", header = F)

# Select "useful" columns
hvs_hits <- hvs_hits[, .(V1, V2, V3)]
colnames(hvs_hits) <- c("read_id", "protein_id", "tax_id")

# Get the virus names
virus_names <- getTaxonomy(ids = hvs_hits[,tax_id], sqlFile = "/scratch/axiom/FAC/FBM/DBC/amalaspi/virome/yarizmen/virome/accessionTaxa.sql", desiredTaxa = c("species"))
# sapply: several requests, then it won't complain about the number of requests
# virus_names <- sapply(hvs_hits[, tax_id], function(x) ncbi_get_taxon_summary(id = x, key = "8f5a7ee70253e037c5641bab48b0d6d74908")$name, USE.NAMES = F, simplify = "vector")
# virus_names <- lapply(hvs_hits[, tax_id], function(x) ncbi_get_taxon_summary(id = x, key = "8f5a7ee70253e037c5641bab48b0d6d74908")$name)
# virus_names <- as.vector(virus_names, mode = "character")
# vapply(hvs_hits[, tax_id], FUN = function(x) ncbi_get_taxon_summary(id = x, key = "8f5a7ee70253e037c5641bab48b0d6d74908")$name, FUN.VALUE = character(length(hvs_hits[, tax_id])))
# virus_names <- ncbi_get_taxon_summary(id = hvs_hits[, tax_id], key = "8f5a7ee70253e037c5641bab48b0d6d74908")$name
# callopts = list(http_version = 0),

# warnings()
# message("Number of taxonomic names found: ", length(virus_names))
hvs_hits[, v_names:=virus_names]

#### Count hits by virus ####

# Count hits per virus
hits_per_virus <- as.data.table(table(hvs_hits[, v_names]))
# Don't lose taxid.
present_virus <- as.data.table(unique(hvs_hits[, tax_id, v_names]))

# Index to match the correct order
index <- match(present_virus[, v_names], hits_per_virus[, V1])
present_virus <- present_virus[, hits:=hits_per_virus[, N[index]]]

# Output the table
write.csv(present_virus, file = output, quote = F, row.names = F)

