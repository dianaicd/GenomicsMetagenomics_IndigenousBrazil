args <- commandArgs(TRUE)
all_pos_file <- args[1]
pos_file <- args[2]
out_name <- args[3]

all_pos <- scan(all_pos_file)
pos <- scan(pos_file)

line_numbers <- match(pos, all_pos)
hom_calls <- c(0, sapply(2:length(line_numbers), function(i) line_numbers[i] - line_numbers[i-1] ) )

write.table(hom_calls, out_name, sep = "", quote = F, row.names = F, col.names = F)