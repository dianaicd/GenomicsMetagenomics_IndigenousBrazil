#!/usr/bin/env Rscript

# September 7 2020.
# Script to save in excel workbooks tables with final viral candidates.

# October 9 2020
# Include effective depth of coverage.

# Actions:
# Create a workbook with a worksheet.

# Running example:
# Rscript --vanilla save2excel_final_candidates_table.R --input sorted_final_candidates_table_all_batches_q0.tsv 
# --quality 0 --output_name ""

#### Libraries ####
library(openxlsx)
library(data.table)
library(optparse)

#### Arguments ####
option_list = list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "Path to the table with the candidates info"),
  make_option(c("-q", "--quality"), type = "character", default = NULL, 
              help = "Mapping quality used"),
  make_option(c("-o", "--output_name"), type = "character", default = "",
              help = "Output name prefix")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

input_table <- opt$input
quality <- opt$quality
output_prefix <- opt$output_name

#### Functions ####
# Function to create the sample data table, it adds the depth of coverage column as well.
create_dt <- function(quality) {
  candidates <- fread(input = input_table)
  
  colnames(candidates) <- c("virus_name", "virus_taxID", "accession_number", "genome_strand_nucleic_acid", 
                            "genome_structure", "genome_length", "sample", "reads_before_rmdup", 
                            "reads_after_rmdup", "average_read_length", "damage_five_prime", "damage_three_prime", 
                            "depth_coverage", "SD_depth_coverage", "breadth_coverage", "effective_depth_coverage", 
                            "SD_effective_depth_coverage")
  
  return(candidates)
}

# Function to write on the excel workbook
create_wb <- function(workbook, dtable, ranges) {
  # add a sheet
  addWorksheet(wb = workbook, sheetName = "Viral candidates")
  # header style
  hs1 <- createStyle(fgFill = "white", halign = "CENTER", textDecoration = "Bold",
                     border = c("top", "bottom", "left", "right"),
                     borderStyle = "thin",
                     fontColour = "black", fontSize = 16)
  
  showGridLines(wb = workbook, sheet = "Viral candidates", showGridLines = FALSE)
  
  setColWidths(wb = workbook, sheet = "Viral candidates", cols = 1:ncol(dtable), widths = "auto")
  # write the data on the sheet
  writeData(wb = workbook, sheet = "Viral candidates", x = dtable, startRow = 1, startCol = 1, 
            headerStyle = hs1, borders = "columns", borderStyle = "thin")
  
  # column style
  # percent
  percent_cols <- which(colnames(dtable) %in% c("breadth_coverage"))
  s <- createStyle(numFmt = "0.00%")
  addStyle(wb = workbook, sheet = "Viral candidates", style = s,  cols = percent_cols, 
           rows = 1:nrow(dtable) + 1, gridExpand = TRUE, stack = T)
  # thousands separator
  comma_cols <- which(colnames(dtable) %in% c("genome_length", "reads_before_rmdup", "reads_after_rmdup"))
  s <- createStyle(numFmt = "COMMA")
  addStyle(wb = workbook, sheet = "Viral candidates", style = s, rows = 1:nrow(dtable) + 1, cols = comma_cols,
           gridExpand = TRUE, stack = T)
  # trim decimal places
  decimal_cols <- which(colnames(dtable) %in% c("average_read_length", "damage_five_prime", 
                                                            "damage_three_prime", "depth_coverage", 
                                                            "SD_depth_coverage", "effective_depth_coverage",
                                                            "SD_effective_depth_coverage")
                        )
  s <- createStyle(numFmt = "#0.00")
  addStyle(wb = workbook, sheet = "Viral candidates", style = s, rows = 1:nrow(dtable) + 1, cols = decimal_cols,
           gridExpand = TRUE, stack = T)
  
  # Borders in selected rows (the ones involved in the same virus)
  for (i in 1:nrow(ranges)){
    style_top <- createStyle(border = "top")
    addStyle(wb = workbook, sheet = "Viral candidates", style = style_top, rows = ranges[i, 1] + 1, cols = 1:ncol(dtable), 
             gridExpand = TRUE, stack = TRUE)
    style_bottom <- createStyle(border = "bottom")
    addStyle(wb = workbook, sheet = "Viral candidates", style = style_bottom, rows = ranges[i, 2] + 1, cols = 1:ncol(dtable), 
             gridExpand = TRUE, stack = TRUE)
  }
  
  # Caption
  captions <- c("virus_name: Virus common name.", "virus_taxID: Virus NCBI taxonomic identifier.", 
  "accession_number: Virus reference genome NCBI accession number.", 
  "genome_strand_nucleic_acid: Nucleic acid that make up the virus genome (DNA or RNA). If the genome is single stranded, it is marked with \"ss-\".", 
  "genome_structure: Type of viral genomic structure (circular or linear genome).", 
  "genome_length: Virus genome length (bp).", "sample: National Museum of Brazil (Museu Nacional, MN) identifier.", 
  "reads_before_rmdup: Number of mapped reads to the specified virus before removing those flagged as duplicates by picardtools MarkDuplicates.",
  "reads_after_rmdup: Number of mapped reads to the specified virus after removing those flagged as duplicates by picardtools MarkDuplicates.",
  "average_read_length: Average length of mapped reads to the specified virus with the specified mapping quality or higher after removing duplicates.",
  "damage_five_prime: Proportion of C to T transitions in the first base of the reads.",
  "damage_three_prime: Proportion of G to A transitions in the last base of the reads.",
  "depth_coverage: Depth of coverage on viral genome (number of unique bases divided by the length of the viral genome).",
  "SD_depth_coverage: Standard deviation of the depth of coverage on the viral genome.",
  "breadth_coverage: Percentage of the viral genome with a depth of coverage >= 1x.",
  "effective_depth_coverage: depth of coverage on positions with coverage >= 1x.",
  "SD_effective_depth_coverage: Standard deviation of the depth of coverage on positions with depth of coverage >= 1x.")

  dt_captions <- data.table(captions)
  # write the data on the sheet
  writeData(wb = workbook, sheet = "Viral candidates", x = dt_captions, startRow = (nrow(dtable) + 2), startCol = 1, 
            colNames = F, borders = "none")

  # format data font type
  s <- createStyle(fontSize = 14, halign = "center", valign = "center")
  addStyle(wb = workbook, sheet = "Viral candidates", style = s, rows = 1:(nrow(dtable) + 1), 
           cols = 1:ncol(dtable), gridExpand = TRUE, stack = T)

  # caption font
  s <- createStyle(fontSize = 14)
  addStyle(wb = workbook, sheet = "Viral candidates", style = s, 
           rows = (nrow(dtable) + 2):(nrow(dtable) + 2 + nrow(dt_captions)), 
           cols = 1:ncol(dtable), gridExpand = TRUE, stack = T)
  
  return(workbook)
}

# Function to select the rows to merge
select_merging_rows <- function(candidates_table){
  # Get the counts to know the legnth of the row ranges
  virus_counts <- table(candidates_table$virus_name)
  
  if (length(virus_counts)==1) {
    # If there is only one type of virus
    row_ranges <- matrix(data = c(1, virus_counts[1]), ncol = 2)
    return(row_ranges)
  } else {
    # Firs range
    row_ranges <- matrix(data = c(1, virus_counts[1]), ncol = 2)
    
    #The rest of ranges
    for (i in c(2:length(virus_counts))) {
      new_range <- c(row_ranges[i - 1, 2] + 1, row_ranges[i - 1, 2] + virus_counts[i])
      row_ranges <- rbind(row_ranges, new_range, deparse.level = 0)
    }
    colnames(row_ranges) <- NULL
    return(row_ranges)
  }
}

#### Create the workbook ####

wb <- createWorkbook()

# Load table
candidates <- create_dt(quality = quality)

# Get row ranges involving same viruses 
row_ranges <- select_merging_rows(candidates_table = candidates)

create_wb(workbook = wb, dtable = candidates, ranges = row_ranges)

# Merge rows involving the same virus
for (i in c(1:6)) {
  for (j in c(1:nrow(row_ranges))) {
    mergeCells(wb, "Viral candidates", cols = i, rows = (row_ranges[j, 1] + 1):(row_ranges[j, 2] + 1) )
  }
}

saveWorkbook(wb, paste0("excel/", output_prefix, "_candidate_virus_q", quality, ".xlsx"), TRUE)