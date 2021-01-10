# Write a dataframe to an Excel sheet

# We need the data, header
# optional: caption, styles


#==============================================================================#
add_new_sheet <- function( df, sheet_name, workbook, column_styles, 
                           caption = F, header_style = F, color_by_sample = F ){
  
  addWorksheet(wb = workbook, sheetName = sheet_name)
  ## Formatting can be applied simply through the write functions 
  showGridLines(wb = workbook, sheet = sheet_name, showGridLines = FALSE)
  setColWidths(wb = workbook, sheet = sheet_name, 
               cols = 2:ncol(df), widths = "auto")
  #--------------
  ## headerStyles
  header_style <- createStyle(fgFill = "white",
                              halign = "CENTER",
                              textDecoration = "Bold",
                              border = c("top", "bottom", "left", "right"),
                              borderStyle = "thin",
                              fontColour = "black", fontSize = 16)
  #----------------
  # Is there a caption?
  # If the caption is embedded in the table, we should remove it as the
  # number format would not be aplied to a column with mix of str and numbers
  # 
  if ( !is.data.frame(caption) ) {
      data_rows <- nrow(df) - ncol(df)
      caption <- df[ data_rows:nrow(df), ]
      df <- df[ 1:data_rows, ]
      na_row <- which( sapply(1:data_rows,
                              function(x) sum(is.na(df[x,]))) == ncol(df) )
      df <- df[ -na_row, ]
      df <- data.frame(apply(df, 2, factor))
  } else {
      data_rows <- nrow(df)
    }
  
  #--------------------------------
  # write table 
  writeData(workbook, sheet = sheet_name, x = df,
            startRow = 1, startCol = 1, headerStyle = header_style,
            borders = "columns", borderStyle = "thin")

  #-------------------
  # add column styles
  add_style_column <- function(column_name, style_list){

    column_index <- which(colnames(df) == column_name)
    style <- do.call(createStyle, style_list)
    
    addStyle(wb = workbook, sheet = sheet_name, 
             style = style,  cols = column_index, 
             rows = 1 : data_rows + 1, 
             gridExpand = TRUE, stack = T) 
    }
  
  lapply(1:length(column_styles),
         function(x) add_style_column( column_name = names( column_styles[x] ),
                                       style_list = column_styles[[x]]
                                       )
         )
  style <- createStyle(fontSize = 14)
  addStyle(wb = workbook, sheet = sheet_name, 
           style = style, 
           rows = 1 : ( data_rows + nrow(caption) )  + 1, 
           cols = 1 : ncol(df),
           gridExpand = TRUE, stack = T)
  #-----------------------------------------------------------------------------#
  # write caption
  # theCaption <- myT1[is.na(myT1$reads_raw),]
  
  writeData(wb = workbook, sheet = sheet_name, 
            x = caption, 
            startRow = nrow(df) + 2,
            startCol = 1,
            colNames = F, 
            borders = "none")
  #----------------------------------------------------------------------------#
  # Optional: Color rows by sample
  if ( color_by_sample ){
    add_color_row <- function( sample ){
      # do_color <- color_rows[ which( samples == sample) ]
      # if( do_color ){
        row_index <- which( df$Sample == sample & df$Library == "All" )
        addStyle(wb = workbook, sheet = sheet_name, 
                 style = style, 
                 rows = row_index + 1, 
                 cols = 1 : ncol(df),
                 gridExpand = TRUE, stack = T)
      # }
    }

    style <- createStyle( fgFill = "#a6e6db" )
    samples <- unique( df$Sample )
    # color_rows <- rep_len( c( F, T ), length.out = length( samples ) )
    
    sapply(samples, function(sample) add_color_row( sample = sample ))
    
  }
}


             