#!/usr/bin/env Rscript

# July 15 2020.
# Script to gather in excel workbooks the output of the Snakefile mystats.

# Actions:
# Create a workbook with a worksheet per sample.
# Tables (MN*_stats.txt) with the output from mystats snakemake are read.
# Depth of coverage is calculated, reference id column is added.

# Running example:
# Rscript --vanilla save2excel.R --refs c(\"AF113323_1\",\"AJ249437_1\",\"AJ717293_1\",
                                          # \"AY083234_1\",\"DQ333427_1\",\"DQ357065_1\",
                                          # \"FN669502_1\",\"HQ340602_1\",\"NC_000883_2\",
                                          # \"NC_001540_1\",\"NC_004295_1\")" --quality 0

#### Libraries ####
library(openxlsx)
library(data.table)
library(optparse)

#### Arguments ####
option_list = list(
  make_option(c("-r", "--refs"), type = "character", default = NULL, 
              help = "reference IDs", metavar = "character"),
  make_option(c("-q", "--quality"), type = "character", default = NULL,
              help = "Mapping quality used")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

reference_ids <- eval(parse(text = (opt$refs)))
quality <- opt$quality

#### Functions ####
# Function to create the sample data table, it adds the depth of coverage column as well.
create_sample_dt <- function(sample_id, quality) {
  sample <- fread(input = paste(sample_id, "/stats/quality_", quality, "/", sample_id, 
                                "_stats.txt", sep = ""))
  
  # CDS_positions <- fread(input = paste0("/scratch/axiom/FAC/FBM/DBC/amalaspi/virome/yarizmen/",
  #                                       "virome/Boto_virome/true_complete/virus_mappings/",
  #                                       "parvovirus/refs_seqs_parvovirus_BMuhlemann/fasta/CDS/",
  #                                       "MN00346/BED_files/refs_parvo.bed"))
  
  sample[, reference_ids := reference_ids]
  
  setnames(sample, 
           names(sample[, 1:8]), 
           c("rds_b4_rmdup", "rds_after_rmdup", "avg_rd_lgth", "BoC", "DoC", "SD_BoC", "SD_DoC", 
           "ref_lgth"),
           skip_absent = T)
  
  # sample[, depth_cov := (rds_after_rmdup*avg_rd_lgth)/ref_lgth]
  
  # sample[, depth_cov_cds := (rds_after_rmdup*avg_rd_lgth)/(ref_lgth - (CDS_positions$V3 - CDS_positions$V2))]
  
  # setcolorder(sample, 
  #             c("reference_ids", "rds_b4_rmdup","rds_after_rmdup","avg_rd_lgth", "genome_cov", 
  #               "ref_lgth", "depth_cov", "dam_3prime", "dam_5prime"))
  
  return(sample)
}

# Function to write on the excel workbook
create_wb <- function(workbook, sample, dtable) {
  # add a sheet
  addWorksheet(wb = workbook, sheetName = sample)
  # header style
  hs1 <- createStyle(fgFill = "white", halign = "CENTER", textDecoration = "Bold",
                     border = c("top", "bottom", "left", "right"),
                     borderStyle = "thin",
                     fontColour = "black", fontSize = 16)
  
  showGridLines(wb = workbook, sheet = sample, showGridLines = FALSE)
  
  setColWidths(wb = workbook, sheet = sample, cols = 1:ncol(dtable), widths = "auto")
  # write the data on the sheet
  writeData(wb = workbook, sheet = sample, x = dtable, startRow = 1, startCol = 1, 
            headerStyle = hs1, borders = "columns", borderStyle = "thin")
  
  # column style
  # percent
  percent_cols <- which(colnames(dtable) %in% c("BoC"))
  s <- createStyle(numFmt = "0.00%")
  addStyle(wb = workbook, sheet = sample, style = s,  cols = percent_cols, 
           rows = 1:nrow(dtable) + 1, gridExpand = TRUE, stack = T)
  # thousands separator
  comma_cols <- which(colnames(dtable) %in% c("rds_b4_rmdup", "rds_after_rmdup", "ref_lgth"))
  s <- createStyle(numFmt = "COMMA")
  addStyle(wb = workbook, sheet = sample, style = s, rows = 1:nrow(dtable) + 1, cols = comma_cols,
           gridExpand = TRUE, stack = T)
  # trim decimal places
  decimal_cols <- which(colnames(sample) %in% c("avg_rd_lgth", "DoC", "SD_BoC", "SD_DoC"))
  s <- createStyle(numFmt = "#0.00")
  addStyle(wb = workbook, sheet = sample, style = s, rows = 1:nrow(dtable) + 1, cols = decimal_cols,
           gridExpand = TRUE, stack = T)
  
  #format data font type
  s <- createStyle(fontSize = 14)
  addStyle(wb = workbook, sheet = sample, style = s, rows = 1:nrow(dtable) + 1, 
           cols = 1:ncol(dtable), gridExpand = TRUE, stack = T)
  
  return(workbook)
}

#### Create the workbook ####

# Get samples' IDs using the dir structure of the snakemake mystats pipeline. 
samples <- list.files(".", include.dirs = T, pattern = "MN")

wb <- createWorkbook()

for (s in samples) {
  create_wb(workbook = wb, 
            sample = s, 
            dtable = create_sample_dt(sample_id = s, quality = quality))
}

saveWorkbook(wb, paste("excel/quality_", quality, ".xlsx", sep = ""), TRUE)