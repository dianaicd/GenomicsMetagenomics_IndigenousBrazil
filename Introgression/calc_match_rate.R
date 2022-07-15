args <- commandArgs(TRUE)
path_tped_obs <- args[1]
path_vars_blocks <- args[2]
path_tped_source <- args[3]
path_ancestral <- args[4]
output_rate <- args[5]


tped_obs <- read.table(path_tped_obs, col.names = c("chr", "id", "x", "pos", "ref", "alt"))
vars_blocks <- read.table(path_vars_blocks, header = T)
tped_obs$ancestral <- toupper( read.table(path_ancestral)$V2 )

tped_obs$derived <- sapply(
    1:nrow(tped_obs),
    function(i)
        ifelse(tped_obs$ref[i] == tped_obs$ancestral[i], tped_obs$alt[i], tped_obs$ref[i])
    )




#path_tped_source <- "Mixe.ascertained.tped"
tped_source <- read.table(path_tped_source, col.names = c("chr", "id", "x", "pos", "ref_source", "alt_source"))
tped_source$id <- NULL

n_block <- 1
chr <- unique(vars_blocks$chr[vars_blocks$n_block == n_block])
pos <- vars_blocks$pos[vars_blocks$n_block == n_block]

require(plyr)

merged_tped <- join(tped_obs, tped_source)
merged_tped <- join(merged_tped[complete.cases(merged_tped),], vars_blocks)


merged_tped$derived_match <- sapply(
    1:nrow(merged_tped),
    function(i)
        merged_tped$derived[i] %in% merged_tped[i, c("ref_source", "alt_source")] 
    )

match_prop <- 
do.call(rbind,
    lapply(
        unique(merged_tped$n_block),
        function(i){
            c(
                i,
                mean(merged_tped$derived_match[merged_tped$n_block == i]),
                unique(merged_tped$prob[merged_tped$n_block == i])   
            )
            }
    )
)


# merged_tped$ref_match <- sapply(
#     1:nrow(merged_tped),
#     function(i)
#         merged_tped$ref[i] %in% merged_tped[i, c("ref_source", "alt_source")] 
#     )

# merged_tped$alt_match <- sapply(
#     1:nrow(merged_tped),
#     function(i)
#         merged_tped$alt[i] %in% merged_tped[i, c("ref_source", "alt_source")] 
#     )

# sum(merged_tped$ref_match | merged_tped$alt_match) / nrow(merged_tped)

# match_prop <- 
# do.call(rbind,
#     lapply(
#         unique(merged_tped$n_block),
#         function(i){
#             c(
#                 i,
#                 mean(merged_tped$ref_match[merged_tped$n_block == i] | merged_tped$alt_match[merged_tped$n_block == i]),
#                 unique(merged_tped$prob[merged_tped$n_block == i])   
#             )
#             }
#     )
# )

write.table(match_prop, output_rate, col.names=F, row.names=F, sep = "\t", quote=F)