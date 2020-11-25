#-------------------------------------------------------------------------------------------------#
# Functions to read and prepare dataframes
read_roh <- function(prefix, sufix,
  file_dir, 
                     chromosomes = seq(1,22),
                     snp_type = "all",
                     change_fid = F,
                     new_id = ''){
  
  roh <- data.frame()
  for (chr in chromosomes) {
    f <- read.table( paste(prefix, chr, sufix, sep = ""), 
                     stringsAsFactors = F,
                     header = T 
                     )
    roh <- rbind(roh, f)  
  }
  
  if ( change_fid ) {
    if( new_id %in% c('1x_fromVictor', '0.5x_fromVictor')){
      for (old_id in unique(roh$IID)){
        index <- match( old_id, fromVictor$old_id)
        roh$IID[ roh$IID == old_id ] <- as.character( fromVictor$new_id[ index ] )
      }
    }else{
      roh$IID <- new_id
    }
  }
  
  roh$FID <- as.character(roh$IID)
  roh$IID <- as.character(roh$IID)
  
  return(roh)
}


prepare_rohs <- function(groups, prefix, sufix, change_fid){
  
  roh_mat <- sapply( 1:length(groups), 
                     function(i) read_roh(prefix = paste(prefix, groups[i], '.chr', sep = ''),
                                          sufix = sufix,
                                          change_fid = change_fid[i], 
                                          new_id = groups[i]))
  
  if(length(groups) > 1){
    roh_list <- apply(roh_mat, 2, function(m) do.call(cbind, m))
    roh_df <- as.data.frame(do.call(rbind, roh_list))
    
  }else{
    roh_df <- as.data.frame( apply(roh_mat, 1, function(m) do.call(cbind, m)) )
  }
  
  for (column in c("KB", "CHR", "POS1", "POS2")) {
    roh_df[, column] <- as.numeric( as.character(roh_df[, column]) )
  }
  
  return(roh_df)
}

# read_roh <- function(file_dir, chromosomes = seq(1,22), snp_type = "all"){
#   roh <- data.frame()
#   for(chr in chromosomes){
#     f <- read.table(paste(file_dir, "chr",chr, ".phased_", snp_type, ".hom", sep = ""), 
#                     header = T)
#     roh <- rbind(roh, f)
#   }
#   return(roh)
# }

get_lengths_counts <- function(roh, breaks_mb, samples){
  
  get_lengths <- function(sample){
    Mb <- roh$KB[ roh$IID == sample ] / 1000
    roh_length <- sapply(1:length(breaks_mb), function(i) sum(Mb[ Mb >= breaks_mb[i] & Mb < breaks_mb[i + 1]]))
    
    data.frame( sample = sample,
                category = breaks_mb,
                length = roh_length)
  }
  
  get_counts <- function(sample){
    Mb <- roh$KB[ roh$IID == sample ] / 1000
    roh_counts <- sapply(1:length(breaks_mb), function(i) length(Mb[ Mb >= breaks_mb[i] & Mb < breaks_mb[i + 1]]))
    
    data.frame( sample = sample,
                category = breaks_mb,
                count = roh_counts)
  }
  roh_counts <- do.call(rbind, lapply(samples, function(sample) get_counts(sample)))
  roh_lengths <- do.call(rbind, lapply(samples, function(sample) get_lengths(sample)))
  lengths_counts <- list( lengths = roh_lengths, counts = roh_counts)
  
  return(lengths_counts)
}


get_sroh_nroh <- function(lengths_counts, samples){
  sroh <- sapply(samples, function(sample) sum(lengths_counts$lengths[, sample]))
  nroh <- sapply(samples, function(sample) sum(lengths_counts$counts[, sample]))
  
  df <- data.frame( sample = samples, 
                    sroh = sroh,
                    nroh = nroh)
  
  return(df)
}

rm_zero_counts_sample <- function(lengths_counts, sample, type, ind_pop){
  sums_with_data <- lengths_counts$lengths[,sample]
  names(sums_with_data) <- seq(1, length(sums_with_data))
  # sums_with_data <- sums_with_data[ sums_with_data != 0 ] 
  
  df = data.frame(length = as.integer(names(sums_with_data)),
                  sum = sums_with_data,
                  type = type,
                  sample = sample,
                  pop = ind_pop$pop[ which(ind_pop$sample == sample) ])
  return(df)
}


merge_roh_lengths <- function(lengths_counts, samples, ind_pop){
  
  build_length_df <- function(lengths_counts, sample, type, ind_pop){
    if(sample %in% colnames(lengths_counts$length)){
      
      lengths_sample <- lengths_counts$lengths[, as.character(sample)]
      names(lengths_sample) <- seq(1, length(lengths_sample))
      
      df = data.frame(length = as.integer(names(lengths_sample)),
                      sum = lengths_sample,
                      type = type,
                      sample = sample,
                      pop = ind_pop$pop[ which(ind_pop$sample == sample) ])
      return(df)
    }
    
  }
  
  length_per_sample <- lapply(samples,
                              function(sample) 
                                do.call(rbind, lapply(names(lengths_counts), 
                                                      function(list_name) 
                                                        build_length_df(lengths_counts = lengths_counts[[list_name]], 
                                                                        sample = sample, 
                                                                        type = list_name, 
                                                                        ind_pop = ind_pop)
                                )
                                )
                              
  )
  
  lengths_merged <- do.call(rbind, length_per_sample)
  
  return(lengths_merged)
}

#-------------------------------------------------------------------------------------------------#
# Functions for plotting
plot_roh_chromosome <- function(roh, chr, header){
  nInd <- length( unique(roh$IID) )
  size <- as.integer(header$size[ header$chr == chr ])
  
  height_chunk <-  1/nInd*9/10
  starting_y <- seq(0, 1 - height_chunk, length.out = nInd)
  ending_y <- seq(height_chunk, 1, length.out = nInd )
  
  # Plot gray background
  plot(bty = "n",
       ylim = c(0, 1), xlim = c(0, 3e8),
       type = "n", x = 0, y = 0, xlab = "Coordinates (Mb)", ylab = NA,
       main = paste('chr', chr),
       axes = F)
  
  rect(xleft = rep(1, nInd), xright = rep(size, nInd),
       ybottom = starting_y, ytop = ending_y, lwd = .2,
       col = "gray90")
  
  # take care of the axes
  individuals <- levels(factor(roh$IID))
  mid_point <- ( ending_y - starting_y )/2
  
  x_ticks <- seq(0, 3e8, 5e7)
  
  axis( 1, at = x_ticks, labels = F )
  text( x = x_ticks,  par("usr")[3], 
        labels = x_ticks/( 1e6 ), pos = 1, xpd = TRUE,
        offset = 1)
  
  text( starting_y + mid_point,  
        labels = individuals, 
        pos = 2, xpd = TRUE, cex = 0.5 )
  
  for ( ind in individuals ){
    
    # Each of the boxes will be colored according to the order of the individuals
    box_num <- as.numeric( as.factor(individuals)[ which(individuals == ind) ] )
    index <- which( roh$IID == ind )
    roh_ind <- roh[index, ]
    
    if(chr %in% roh_ind$CHR){
      # Add ROH
      rect(xleft = roh_ind$POS1[ roh_ind$CHR == chr ],
           xright = roh_ind$POS2[ roh_ind$CHR == chr ],
           ybottom = starting_y[ box_num ],
           ytop = ending_y[ box_num ], lwd = .2,
           col = "darkcyan") 
    }
  }
}

plot_roh_ind <- function(roh, ind, plot_dir, header){
  
  # roh <- roh[roh$IID == ind, ]
  header$chr <- sapply(header$V2, function(x) strsplit(x, ":")[[1]][2])
  header$size <- sapply(header$V3, function(x) strsplit(x, ":")[[1]][2])
  
  fig_name <- paste(plot_dir, "ROH_",ind,".png", sep = "")
  png(fig_name,
      height = 6, width = 12, res = 250, units = "in")
  layout(matrix(seq(1, 24), byrow = T, nrow = 4, ncol = 6),
         widths = rep(2,6), heights = rep(1,4))
  
  for(chr in 1:22){
    
    
  }
  dev.off()
}

plot_roh_all_ind <- function(roh, plot_dir, header, ncols = 6, nrows = 4, width_box = 2, height_box = 2){
  
  # roh <- roh[roh$IID == ind, ]
  header$chr <- sapply(header$V2, function(x) strsplit(x, ":")[[1]][2])
  header$size <- sapply(header$V3, function(x) strsplit(x, ":")[[1]][2])
  
  fig_name <- paste(plot_dir, "ROH_all_ind.png", sep = "")
  png(fig_name,
      height = 2*nrows, width = 3*ncols, res = 350, units = "in")
  layout(matrix(seq(1, ncols*nrows), byrow = T, nrow = nrows, ncol = ncols),
         widths = rep(width_box, ncols), heights = rep(height_box, nrows))
  
  for(chr in 1:22){
    plot_roh_chromosome(roh = roh, chr = chr, header = header)
    
  }
  dev.off()
}


plot_roh_lengths_sum <- function(lengths_counts, my_legend, 
                                 breaks_mb,
                                 xlim = c(0, 8), x_tick = seq(0,8),
                                 ylim = c(0, 125),
                                 add_legend = F,
                                 legend_coords = c(5, 100),
                                 main = 'ROH per genome', cex.main = 1.25
){
  
  plot(1,1, xlim = xlim, ylim = ylim, axes = F, 
       main = main, cex.main = cex.main,
       xlab = 'ROH length category (Mb)',
       ylab = 'Sum of lengths (Mb)', type = 'n')
  
  axis(1, at = x_tick)
  axis(2)
  grid()
  
  for( sample in rev( my_legend$samples ) ){
    
    my_sroh <- lengths_counts$lengths[, sample]
    # my_sroh <- my_sroh[ my_sroh != 0]
    index <- which(my_legend$samples == sample)
    
    points(x = breaks_mb[which(my_sroh != 0) + 1], y = my_sroh[ my_sroh != 0], 
           type = 'b',
           col = as.character(my_legend$color[ index ]),
           pch = my_legend$pch[ index ],
           lwd = 2)
    
  }
  
  if( add_legend ){
    legend(legend_coords[1], legend_coords[2], legend = my_legend$legend, 
           pch = my_legend$pch, 
           col = as.character(my_legend$color),
           bty = 'n',
           pt.lwd = 2)
    
  }
}

plot_roh_lengths_counts <- function(lengths_counts, my_legend, 
                                 breaks_mb,
                                 xlim = c(0, 8), x_tick = seq(0,8),
                                 ylim = c(0, 125),
                                 add_legend = F,
                                 legend_coords = c(5, 100),
                                 main = 'ROH per genome', cex.main = 1.25
){
  
  plot(1,1, xlim = xlim, ylim = ylim, axes = F, 
       main = main, cex.main = cex.main,
       xlab = 'ROH length category (Mb)',
       ylab = 'Number of ROH', type = 'n')
  
  axis(1, at = x_tick)
  axis(2)
  grid()
  
  for( sample in rev( my_legend$samples ) ){
    
    my_nroh <- sapply(lengths_counts$counts, function(l) unlist(l)[sample])
    # my_sroh <- my_sroh[ my_sroh != 0]
    index <- which(my_legend$samples == sample)
    
    points(x = breaks_mb[which(!is.na(my_nroh)) + 1], y = my_nroh[!is.na(my_nroh)], 
           type = 'b',
           col = as.character(my_legend$color[ index ]),
           pch = my_legend$pch[ index ],
           lwd = 2)
    
  }
  
  if( add_legend ){
    legend(legend_coords[1], legend_coords[2], legend = my_legend$legend, 
           pch = my_legend$pch, 
           col = as.character(my_legend$color),
           bty = 'n',
           pt.lwd = 2)
    
  }
}

plot_nroh_sroh <- function( groups, prefix, sufix, change_fid,
                            breaks, ind_pop, colors, add_label,
                            xlim = c(0, 900),
                            ylim = c(0, 900)){
  
  
  roh_all <- prepare_rohs(groups = groups, 
                          prefix = prefix,
                          sufix = sufix,
                          change_fid = change_fid)
  
  samples <- unique(roh_all$IID)
  
  lengths_counts <- get_lengths_counts(roh = roh_all, 
                                       breaks = breaks, 
                                       samples = samples)
  
  sroh_nroh <- get_sroh_nroh( lengths_counts = lengths_counts,
                              samples = samples )
  sroh_nroh <- join(sroh_nroh, ind_pop, by = 'sample')
  
  
  
  roh_plot <- ggplot(sroh_nroh, aes(x = sroh, y = nroh, color = pop)) +
    geom_point(size = 4) +
    coord_cartesian( xlim = xlim,
                     ylim = ylim) +
    labs(x = 'Sum of ROH lengths in the genome (Mb)',
         y = 'Number of ROHs per genome') +
    scale_color_manual(values = colors)
  
  if( add_label ){
    roh_plot <- roh_plot + geom_text( label = sroh_nroh$sample )
  }
  
  return( roh_plot )
  
}

gather_plot <- function(){
  
  sufix <- '.phased_all_roh.hom'
  all_plot <- plot_nroh_sroh( groups,
                              prefix = prefix,
                              sufix = sufix,
                              change_fid,
                              breaks = seq( 0, 15),
                              ind_pop = ind_pop, 
                              colors = colors, add_label = F,
                              xlim = c(0, 900),
                              ylim = c(0, 900))
  
  all_plot_label <- plot_nroh_sroh( groups,
                                    prefix = prefix, sufix = sufix,
                                    change_fid,
                                    breaks = seq( 0, 15),
                                    ind_pop = ind_pop, 
                                    colors = colors, add_label = T,
                                    xlim = c(0, 900),
                                    ylim = c(0, 900))
  
  
  rmTrans_plot <- plot_nroh_sroh( groups = groups,
                                  prefix = prefix, sufix = sufix,
                                  change_fid = change_fid,
                                  breaks = seq( 0, 35),
                                  ind_pop = ind_pop, 
                                  colors = colors, add_label = F,
                                  xlim = c(0, 900),
                                  ylim = c(0, 900))
  
  rmTrans_plot_label <- plot_nroh_sroh( groups = groups,
                                        prefix = prefix, sufix = sufix,
                                        change_fid = change_fid,
                                        breaks = seq( 0, 35),
                                        ind_pop = ind_pop, 
                                        colors = colors, add_label = T,
                                        xlim = c(0, 900),
                                        ylim = c(0, 900))
  
  leyenda <- get_legend(all_plot)
  
  all_plot <- all_plot + theme(legend.position = 'none')
  rmTrans_plot <- rmTrans_plot + theme(legend.position = 'none')
  all_plot_label <- all_plot_label + theme(legend.position = 'none')
  rmTrans_plot_label <- rmTrans_plot_label + theme(legend.position = 'none')
  
  plot_grid(all_plot, rmTrans_plot,  leyenda,
            all_plot_label, rmTrans_plot_label,
            ncol = 3, rel_widths = c(3,3,1))
}


plot_lengths_pop_sample <- function(df, low_color = "#ffab00", high_color = "#008aa1",
                                    ylim = c(0, 300), ncol = 2, 
                                    colors = c("#ff616b", "#ff616b", "#328e13", "#172713", "#a6e6db")){
  
  ggplot(df, aes(x = sample, y = sum) ) +
    geom_bar(stat = 'identity',
             position =  position_dodge2(preserve = "single"),
             aes(fill = length)) +
    labs(x = '',
         y= 'Sum (Mb)') +
    coord_cartesian(ylim = ylim) +
    facet_wrap(. ~ pop, 
               scales = 'free_x',
               ncol = ncol) +
    scale_fill_gradientn(#low = low_color,
                        #high = high_color,
                        colours = colors,
                        aesthetics = "fill", 
                        space = "Lab",
                        name = 'Length (Mb)') +
    theme(axis.text.x = element_text(angle = 90, size = 10))
}


points_length_sum <- function(lengths_counts, sample, color ){
  
  sums_with_data <- lengths_counts$lengths[,sample]
  names(sums_with_data) <- seq(1, length(sums_with_data))
  sums_with_data <- sums_with_data[ sums_with_data != 0 ] 
  points(x = names(sums_with_data),
         y = sums_with_data,
         type = 'b', col = color, lwd = 2)
  
}



plot_sroh_nroh_trajectory <- function(all_sroh_nroh, xlim, ylim, myColors){
  
  populations <- unique(all_sroh_nroh$pop)
  if( ! 'pch' %in% colnames(all_sroh_nroh) ){
    all_sroh_nroh$pch <- as.numeric(factor(all_sroh_nroh$type))
  }
  
  layout(matrix(seq(1, nrow*ncol), nrow = nrow, byrow = T))
  
  for(pop in populations){ 
    pop_sroh_nroh <- all_sroh_nroh[all_sroh_nroh$pop == pop,]
    pop_sroh_nroh$color <- as.character(myColors[as.numeric(factor(pop_sroh_nroh$sample))])
    plot(1, 1, type = 'n', xlim = xlim, 
         ylim = ylim, main = pop,
         axes = F, xlab = 'Sum of ROH lengths per genome (Mb)',
         ylab = 'Number of ROH per genome')
    axis(1)
    axis(2)
    grid()
    samples <- unique(pop_sroh_nroh$sample)
    
    for(sample in samples){
      
      index <- which(pop_sroh_nroh$sample == sample)
      
      points(pop_sroh_nroh$sroh[index], 
             pop_sroh_nroh$nroh[index], 
             type = 'b', pch = pop_sroh_nroh$pch[index],
             col = pop_sroh_nroh$color[index],
             lwd = 2)
    }
  }
  
  myLegend <- unique(pop_sroh_nroh[, c('type', 'pch')])
  plot(1:10, 1:10, type = 'n', main = NA, xlab = NA, ylab = NA, axes = F)
  legend(2, 8,legend = myLegend$type, pch = myLegend$pch,
         bty = 'n')
}
#-------------------------------------------------------------------------------------------------#