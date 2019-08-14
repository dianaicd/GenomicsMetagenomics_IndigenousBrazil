library(contamMix)

args <- commandArgs(trailingOnly = TRUE)
load(args[1])
print.contamMix(res)