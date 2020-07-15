###############################################################################
# Get clean tables with values and data types you can still use in R 

# Function to scale values from 0 to 100 to 0 to 1
percent_to_fraction <- function( myTable ){
  myTable$trim_prop       <- myTable$trim_prop / 100
  myTable$duplicates_prop <- myTable$duplicates_prop / 100
  myTable$endo_prop       <- myTable$endo_prop / 100
  myTable$endo_final_prop <- myTable$endo_final_prop / 100
  return(myTable)
}

# Samples' order by deppth of coverage
get_samples_order <- function( myTable ){
  sample_DoC <- unique(myTable[, c("Sample", "DoC_endogenous")])
  samples_order <- unique(myTable$Sample[ order(myTable$DoC_endogenous) ])
  return( samples_order )
}

order_samples <- function( myTable ){
  samples_order <- get_samples_order( myTable = myTable )
  myTable$Sample <- factor(myTable$Sample, 
                           levels = samples_order,
                           ordered = T)
  myTable <- myTable[order(myTable$Sample, decreasing = T),]
  return( myTable )
}
# Reduce and order table to columns listed in caption
reduce_order_columns <- function( myTable, myColnames ){
  # samples_order <- get_samples_order( myTable = myTable )
  column_order <- myColNames$new_name[ myColNames$new_name %in% colnames(myTable) ]
  
  # myTable$Sample <- factor(myTable$Sample, 
  #                       levels = samples_order,
  #                       ordered = T)
  # myTable <- myTable[order(myTable$Sample, decreasing = T),]
  myTable <- myTable[, column_order]
  myTable <- myTable[, colnames(myTable) %in% myColNames$new_name ]
  return( myTable )
}
# Translate table names to those you decided for the caption
change_name <- function(name){
  new_name <- myColNames$new_name[ which( myColNames$orig_name == name) ]
  ifelse( length(new_name) == 0, return(name), return(new_name) )
}
#-----------------------------------------------------------------------------#
# Merge stats from depth and stats files into supplementary table 1 
merge_and_clean_stats <- function( st, dp ){
  myStats <- join( st, dp, by = "SM" )
  myStats <- myStats[ order( myStats$endo_final_prop ), ]
  myStats$SM <- factor( myStats$SM, levels = unique( myStats$SM ), ordered = T )
  
  message( "Please make sure that the paths to 
          split samples MN0008, MN01701 and MN1943 are up to date." )
  myStats <- handle_sample_splits( myStats )
  #---------------------------------------------------------------------------#
  # go from 0 - 100 to 0 - 1 for proportions and further formatting
  myStats <- percent_to_fraction( myStats )
  # format the 95CI of Ry
  myStats$low_ry <-  sapply( strsplit( myStats$X95.CI, "-" ), 
                            function(x) as.numeric( x[1] ) )
  myStats$high_ry <-  sapply( strsplit( myStats$X95.CI, "-" ), 
                             function(x) as.numeric( x[2] ) )
  
  return(myStats)
}
#------------------------------------------------------------------------------#
# Function to do take care of the mess for splitting the sample into  a mix of 
# libraries or separate libraries
split_sample_groups <- function( sample,
                                 libs1, libs2,
                                 new_sm1, new_sm2, 
                                 new_lb1, new_lb2,
                                 path_depth_lib, lib_st ){
  
  myDepth <- read.csv(path_depth_lib)
  
  myDepth$Sample <- sub(".depth.txt", "", myDepth$Sample)
  
  subset_stats <- function(lib_st, sample, libs, new_sample, lib_name){
    sample_stats <- lib_st[lib_st$SM == sample,]
    
    subset_index <-  sample_stats$LB %in% libs
    mySubset_stats <- data.frame(
      SM            = paste(sample, lib_name, sep = "_"),
      SM_new        = new_sample,
      reads_raw     = sum(sample_stats$reads_raw[subset_index]),
      reads_trim    = sum(sample_stats$reads_trim[subset_index]),
      mapping       = sum(sample_stats$mapping[subset_index]),
      duplicates    = sum(sample_stats$duplicates[subset_index]),
      mapping_final = sum(sample_stats$mapping_final[subset_index])
    )
    mySubset_stats$trim_prop       <- mySubset_stats$reads_trim/mySubset_stats$reads_raw*100
    mySubset_stats$endo_prop       <- mySubset_stats$mapping/mySubset_stats$reads_raw*100
    mySubset_stats$duplicates_prop <- mySubset_stats$duplicates/mySubset_stats$mapping*100
    mySubset_stats$endo_final_prop <- mySubset_stats$mapping_final/mySubset_stats$reads_raw*100
    
    
    return(mySubset_stats)
  }
  
  set1 <- subset_stats(lib_st = lib_st, sample = sample, libs = libs1,
                       new_sample = new_sm1, lib_name = new_lb1)
  set2 <- subset_stats(lib_st = lib_st, sample = sample, libs = libs2,
                       new_sample = new_sm2, lib_name = new_lb2)
  
  sample_stats <- rbind(set1, set2)
  myDepth$SM <- myDepth$Sample
  
  mySample <- join(sample_stats, myDepth, by = "SM")
  mySample$SM <- mySample$SM_new
  mySample$Sample <- mySample$SM_new
  mySample$SM_new <- NULL
  mySample$LB <- NULL
  mySample$Library <- NULL
  return(mySample)
}
#------------------------------------------------------------------------------#
# Special split of results for MN0008 (MN0008_non_U, MN0008_L3U),
# MN01701 (remove mtCapture), MN1943 (remove mtCapture) and God knows what else

# Add MN0008_L3U and MN0008_non_U,
# Add MN1943_L1 and MN01701_L1L2 (remove mtCapture)
handle_sample_splits <- function( myStats ){
  lib_st <- read.csv(paste("~/Projects/Botocudos/Files/Summaries/2019_11_21/",
                           data_type, "/library_stats.hg19.csv", sep = ""))
  mn8 <- split_sample_groups(sample = "MN0008",
                             libs1 = c("L3U"), 
                             libs2 = c("L1", "L2"),
                             new_sm1 = "MN0008_L3U",
                             new_sm2 = "MN0008_non_U", 
                             new_lb1 = "L3U",
                             new_lb2 = "non_U",
                             path_depth_lib = paste("~/Projects/Botocudos/Files/Summaries/2019_11_21/",
                                                    data_type, "/MN0008_bylib.depth.txt", sep = ""),
                             lib_st = lib_st)
  
  mn1943 <- split_sample_groups(sample = "MN1943",
                                libs1 = c("L1"), libs2 = c("mtCapture"),
                                new_sm1 = "MN1943",
                                new_sm2 = "MN1943_mtCapture", 
                                new_lb1 = "L1",
                                new_lb2 = "mtCapture",
                                path_depth_lib = paste("~/Projects/Botocudos/Files/Summaries/2019_11_21/",
                                                       data_type,"/MN1943_bylib.depth.txt", sep = ""),
                                lib_st = lib_st)
  
  mn01701 <- split_sample_groups(sample = "MN01701", 
                                 libs1 = c("L1", "L2"),
                                 libs2 = c("mtCapture"),
                                 new_sm1 = "MN01701",
                                 new_sm2 = "MN01701_mtCapture", 
                                 new_lb1 = "L1L2",
                                 new_lb2 = "mtCapture",
                                 path_depth_lib = paste("~/Projects/Botocudos/Files/Summaries/2019_11_21/",
                                                        data_type, "/MN01701_bylib.depth.txt", sep = ""),
                                 lib_st = lib_st)
  
  
  myStats <- myStats[!myStats$SM == "MN1943", ]
  myStats <- myStats[!myStats$SM == "MN01701", ]
  myStats <- rbind(myStats, mn8)
  myStats <- rbind(myStats, mn1943[mn1943$SM == "MN1943",])
  myStats <- rbind(myStats, mn01701[mn01701$SM == "MN01701",])
  
  return(myStats)
}

#------------------------------------------------------------------------------#
# Merge stats from depth and stats files into supplementary table 2 (libraries)

merge_stats_libs <- function( sm_stats, lib_stats, myColNames, backbone, 
                              bamdamage_path ){
  
  # Merge the previously obtained statistics per sample,
  # which would be the values indicated for the libraries "All"
  # Make empty columns if needed to concatenate the tables
  sm_stats$LB  <-  "All"
  lib_stats <- percent_to_fraction(lib_stats)
  missing_cols <- colnames(sm_stats)[ ! colnames( sm_stats ) %in% colnames( lib_stats ) ]
  lib_stats[, missing_cols] <- sapply(missing_cols,
                                      function(x) x = numeric(nrow(lib_stats)))
  
  assign_coverage <- function(sample){
    myCoverage <- sm_stats$AvgReadDepth_tot[sm_stats$Sample == sample]
    return(myCoverage)
  }
  # We need the depth of coverage to order the samples
  lib_stats$AvgReadDepth_tot <- sapply(lib_stats$SM, 
                                       function(x) assign_coverage(x))
  myT2 <- rbind( sm_stats, lib_stats )
  myT2$Sample <- myT2$SM
  myT2$AvgReadLength <- sapply(1:nrow(myT2), 
                               function(x) mean_length(sample = myT2$Sample[x],
                                                       library = myT2$LB[x],
                                                       bamdamage_path = bamdamage_path))
  
  myT2 <- stats_to_readable_tables(myStats = myT2, backbone = backbone,
                                   myColNames = myColNames )
  
  message("Please verify the factor levels for the Library column: mtCapture, L3U, L2, L1, All")
  myT2$Library <- factor( myT2$Library, levels = c( "mtCapture", "L3U", 
                                                    "L2", "L1", "All" ),
                          ordered = T )
  
  myT2    <- myT2[ order( myT2$Sample, myT2$Library, decreasing = T ), ]
  myT2$DoC_endogenous <- NULL
  
  # To do: add read length and damage
  message(paste("Please double check path to bamdamage files: ", 
                bamdamage_path))
  
  return(myT2)
}
################################################################################
# Format Supplementary table 1 so that it is more "human-readable"
# Merge metadata (backbone) to table (myT1)
# Append to the table the description of the columns
stats_to_readable_tables <- function(myStats, 
                                     backbone, 
                                     myColNames, 
                                    add_caption = F){
  
  myTable <- myStats
  

  
  # Ry values are reported as text with their confidence intervals
  # Ry is in ST1 but not in ST2
  if ( "R_y" %in% colnames( myTable ) ){
    myTable$R_y <-  paste(round(myTable$R_y, 3), " (", round(myTable$low_ry, 3),
                          "-",  round(myTable$high_ry, 3),
                          ")", sep = "")
    myTable$low_ry <- NULL
    myTable$high_ry <- NULL
  }
  
  # Add metadata, change column names, 
  # drop columns not present in the caption and order rows and columns
  myTable <- join(backbone, myTable, by = "Sample")
  
  colnames(myTable) <- sapply(colnames(myTable), 
                           function(name) change_name(name) )
  
  myTable <- reduce_order_columns( myTable = myTable,
                                   myColnames = myColnames )
  myTable <- order_samples( myTable = myTable )
  
  # Append empty row and caption if requested
  if(add_caption){
    NA_row <- sapply(myTable[1,], function(x) return(NA))
    myCaption <- get_caption( myColNames = myColNames,
                              myTable = myTable)
    
    myTable <- rbind(myTable, NA_row)
    myTable <- rbind(myTable, myCaption)
  }
  return(myTable)
}

get_caption <- function(myColNames, myT1){
  
  myColNames <- myColNames[ myColNames$new_name %in% colnames(myT1), ]
  
  myCaption <- as.data.frame( sapply( 1:ncol(myT1), 
                                      function(x) 
                                        colnames(myT1)[x] = character(ncol(myT1))
  )
  )
  colnames(myCaption) <- colnames(myT1)
  myCaption[,1] <- unique( paste(myColNames$new_name, 
                                 myColNames$Description, sep = ": ") )
  return(myCaption)
}

################################################################################
# Gather length and damage statistics
read_length_file <- function( sample, bamdamage_path, library = F ){
  if ( library == "All" ){
    library <- F
  }
  if(is.character(library)){
    path_length <- paste(bamdamage_path, "/", sample, "/", 
                         library, "/library_bamdamage/", library,
                         ".hg19.length.csv", sep = "")
  }else{
    path_length <- paste(bamdamage_path, "/", 
                         sample, "/", sample, ".hg19.length.csv", sep = "")  
  }
  myLength <- read.csv(path_length)
  return( myLength )
}

read_damage_file <- function( sample, bamdamage_path, library = F ) {
  if ( library == "All" ){
    library <- F
  }
  if(is.character(library)){
    path_dam <- paste(bamdamage_path, "/", sample, "/", 
                         library, "/library_bamdamage/", library,
                         "..hg19.dam_5prime.csv", sep = "")
  }else{
    path_dam <- paste(bamdamage_path, "/", 
                         sample, "/", sample, ".hg19.dam_5prime.csv", sep = "")  
  }

  damage_5 <- read.csv(path_dam)
  path_dam <- sub("5prime", "3prime", path_dam)
  damage_3 <- read.csv(path_dam)
  
  damage <- list(five_prime = damage_5, three_prime = damage_3)
  return( damage )
}

mean_length <- function(sample, bamdamage_path, library = F){
  myLength <- read_length_file( sample, bamdamage_path, library )
  avg_length <- sum( myLength$counts * myLength$length ) / sum( myLength$counts )
  return(avg_length)
}

damage_first_base <- function(  sample, bamdamage_path, library,
                                end = 5, orig_base = "C", obs_base = "T",
                                pos_base = 1 ){
  
  damage <- read_damage_file( sample, bamdamage_path, library )
  mutation <- paste(orig_base, obs_base, sep = "..")
  if( end == 3 ){
    freq_mutation <- damage$five_prime[ pos_base, mutation ]
  }else if( end == 5 ){
    freq_mutation <- damage$three_prime[ pos_base, mutation ]
  }
  return(freq_mutation)
}

################################################################################
# Clean and make contamination tables

#------------------------------------------------------------------------------#
# ContaminationX

# Start out with a dataframe with "Sample", "Library" and "Sex" columns
read_contaminationX_file <- function(contaminationX_path, sample, library, 
                                     min_depth = 2){
  empty_df <- function(error){
    # print(error)
    myDF <- data.frame(method = c("One-cns", "Two-cns"),
                       estimate = c(NA, NA), lb = c(NA, NA), ub = c(NA, NA),
                       err = c(NA, NA), nSites = c(NA, NA))
    return(myDF)
  }
  result_path <- paste(contaminationX_path, sample,".hg19_", library, 
                       "_md", min_depth, ".result", sep = "")
  estimates_X <- tryCatch(estimates_X <- read.table(result_path),
                          error =  function(error) empty_df(error),
                          finally = function() return(estimates_X))
  colnames(estimates_X) <- c("method", "estimate", "lb", "ub", "err", "nSites")
  estimates_X$Sample <- sample
  estimates_X$Library <- library
  estimates_X$depth <- min_depth
  return(estimates_X)
}

prepare_contaminationX_table <- function(contaminationX_path, samples_libs,
                                         samples_sex){
  message(paste( "Please verify this path is up to date:", contaminationX_path ))
  samples_estimates <- lapply(1:nrow(df), 
                              function(i) 
                                read_contaminationX_file( contaminationX_path = contaminationX_path,
                                                          sample = samples_libs$Sample[i],
                                                          library = samples_libs$Library[i],
                                                          min_depth = 3
                                )
  )
  
  samples_estimates <- do.call(rbind, samples_estimates)
  
  samples_estimates <- join(samples_estimates,
                            samples_sex[, c("Sex", "Sample")], 
                            by = "Sample")
  males_estimate <- samples_estimates[ ! is.na(samples_estimates$estimate) 
                                       & samples_estimates$Sex 
                                       %in% c( "XY", "Not Assigned", 
                                               "consistent with XY but not XX")
                                       & samples_estimates$method == "Two-cns", 
                                       ]
  
  males_estimate$estimate <- paste(percent(males_estimate$estimate),
                                   " (",
                                   percent(males_estimate$lb), 
                                   " - ",
                                   percent(males_estimate$ub),
                                   ")", sep = "")
  
  colnames(males_estimate) <- c("X_method", "contaminationX_estimate", 
                                "X_low_est", "X_upper_est",
                                "X_error", "contaminationX_num_sites",
                                "Sample", "Library", 
                                "X_min_depth",
                                "Sex")
  
  males_estimate <- males_estimate[, c("Sample",
                                       "Library",
                                       "contaminationX_estimate", 
                                       "contaminationX_num_sites",
                                       "X_min_depth"
                                      )
                                  ]
  return(males_estimate)
}

#------------------------------------------------------------------------------#
# schmutzi
get_schmutzi_table <- function(schmutzi_path){
  message(paste("Please verify that the path is up to date: ", schmutzi_path))
  schmutzi_files <- list.files(schmutzi_path, pattern =".est")
  schmutzi <- data.frame()
  
  for(f in schmutzi_files){
    tmp_schmutzi <- read.table(paste(schmutzi_path, f, sep =""))
    colnames(tmp_schmutzi) <- c("schmutzi_estimate", "low", "up")
    tmp_schmutzi$Sample <- sub(".final.cont.est", "", f)
    schmutzi <- rbind(schmutzi, tmp_schmutzi)
  }
  schmutzi$schmutzi_estimate <- paste(percent(schmutzi$schmutzi_estimate),
                             " (",
                             percent(schmutzi$low),
                             " - ",
                             percent(schmutzi$up),
                             ")",
                             sep = "")
  schmutzi[, c("low", "up")] <- NULL
  return(schmutzi)
}

#------------------------------------------------------------------------------#
# contamMix
# We can either 
# load a table with the MAP and also annotate the number of reads used, or
# load the output from the chains and extract a fraction of it (to plot distributions)
prefix <- "~/Projects/Botocudos/Files/Contamination/2020_01_24/contamMix/"

rmTrans <- c("all", "rmTrans")

# We ran the chains for 100k iterations, 
# we will subsample a fraction of these samples to plot the posterior distribution

subsample_contammmix <- function(folder_path, sample, rmTrans, 
                                 size = 100, lib = "All",
                                 burnIn = 0.1){
  
  myRda <- paste( folder_path, sample, "/", sample,
                  ".hg19_", lib, "_", rmTrans, ".Rdata", sep = "")
  
  print(myRda)
  load(myRda)
  notBurnIn <- as.integer(burnIn*size)
  #remove burnin
  tiny <- c(res$chains[[1]][size - ((size- notBurnIn):0),2])
  
  result <- data.frame(estimate = tiny, damage = rmTrans, lib = lib,
                       sample = sample)
  return(result)
}

# contammix <- data.frame()
# for(s in samples){
#   libs <- (myT2$LB[myT2$SM == s])
#   for(l in libs){
#     for(r in rmTrans){
#       tmp_contam <- subsample_contammmix(prefix, s, r, 100000, l)
#       contammix <- rbind(contammix, tmp_contam)
#     }
#   }
# }


#------------------------------------------------------------------------------#
# Order libraries and save (or load) results


get_contamMix_table <- function(contamMix_path,
                                map_path, 
                                samples, 
                                rmTrans = c("all", "rmTrans")){
  
  get_annotation <- function(contamMix_path, sample, library, rmTrans){
    myRda <- paste(contamMix_path, sample, "/", sample, ".hg19_", library, "_",
                   rmTrans, ".Rdata", sep = "")
    load(myRda)
    nReads <- dim(res$mnMatrix)[1]
    
    annot <- data.frame(num_reads_contamMix = nReads,
                        damage = as.character(rmTrans), 
                        Library = as.character(library),
                        Sample = as.character(sample))
    return(annot)
  }
  
  
  annotations <- lapply(samples, function(sample) 
    do.call(rbind, 
            lapply(rmTrans, function(mutation_type)
              get_annotation(contamMix_path = contamMix_path,
                             sample = sample, 
                             library = "All", 
                             rmTrans = mutation_type)
            )
    )
    
  )
  
  annot_mito <- do.call(rbind, annotations)
  
  MAP <- read.csv(map_path, header = F)
  colnames(MAP) <- c("Sample", "Library", "damage","map", "low", "high")
  MAP$Library <- sub("S1", "L1", MAP$Library)
  
  contamMix_table <-  join(annot_mito, MAP, by = c("Sample", "Library", "damage"))
  contamMix_table$contamMix_estimate <- paste(percent(1 - contamMix_table$map,
                                                      accuracy = 0.1),
                                              " (",
                                              percent(1 - contamMix_table$high,
                                                      accuracy = 0.1),
                                              " - ",
                                              percent(1 - contamMix_table$low,
                                                      accuracy = 0.1),
                                              ")",
                                              sep = "")
  
  contamMix_table[,c("map", "low", "high")] <- NULL
  return(contamMix_table)
}